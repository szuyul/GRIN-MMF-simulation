function [ T, NMode, pq_map, betas, E, refractive_idx ]...
              = GRIN_MMF_simTM_HG( lambda, D, NA, Length, Rho, Theta, img_size, N )
%GRIN_MMF_simTM_HG simulates the scalar HG modes, a set of Hermite Gaussian functions, 
% of a GRIN MMF given the specifications and geometry. The HG solutions approximate well to propagation-invariant modes (PIMs.)
% The code is the implementation of theoretical GRIN MMF model:
% Mahdieh B. Shemirani et al., "Principal Modes in Graded-Index Multimode Fiber in Presence of Spatial- and Polarization-Mode Coupling,"
%   J. Lightwave Technol. 27(10), 2009.
%
% The fiber transmission of a deformed MMF is simulated based on a perturbation theory by calculating a coupling matrix, 
% while the HG mode profiles do not change.
%
% [ T, pq_map, betas, E, img_size, refractive_idx ]...
%            = GRIN_MMF_simTM_HG( lambda, D, Length, Rho, Theta, N, n_count_modes ) 
% 
% outputs:
% T is the transmission matrix of a straight or deformed GRIN MMF in HG representation at a certain wavelength
% NMode is the number of counted modes for a GRIN MMF
% pq_map is a 3 by n_count_modes matrix, with each column as the characteristics of one particular mode
%                                      1st row as the angular momentum, p + q + 1
%                                      2nd row as the p indices (# of nulls in x)
%                                      3rd row as the q indices (# of nulls in y)
% beta is a 1 by n_modes matrix, the propagation constants (unit: m) of each HG mode
% E is a N by N by n_modes matrix, the transverse electric scalar field distribution of each HG mode in one polarization
% refractive_idx is the refractive index profile over the img_size
%
% inputs:
% lambda is the wavelength (unit: m)
% D is the fiber core diameter of the MMF (unit: m)
% Length is an array of the segmental length of the MMF (unit: m)
% Rho is an array of the bending radius (unit: m) of each defined segment. If the segment is straight, Rho = inf
% Theta is an array of the orientation of the bending of the MMF (unit: rad.)
% img_size is the physical size (unit: m) of the transverse field images
% N is the image dimension, e.g., 32
%
%
% 2020 Szu-Yu Lee
% BLCTO at Nokia Bell Labs

%% GRIN MMF parameter initialization (unit in m)
k = 2*pi/lambda;
a = D/2; 
n_0 = 1.4875; 
% https://www.fiberoptics4sale.com/blogs/archive-posts/95048070-basic-optics-for-optical-fiber
% n_clad = 1.47; 
n_clad = sqrt(n_0^2 - NA^2);
Dn = (n_0 - n_clad)/n_0;
beta_low = n_0*k*cos(asin(NA/n_0)); % lower bound on beta considering max acceptance angle in MMF
wo = sqrt(D/(k*n_0*sqrt(2*Dn))); % wo is the waist of the fundamental HG mode (unit: m)
M = a*k*n_0*sqrt(Dn/2);

range = linspace(-img_size/2, img_size/2, N);
refractive_idx = n_0*sqrt(1 - 2*Dn*(range/a).^2);
[x, y] = meshgrid(range);
E0 = exp(-(x.^2+y.^2)/wo^2); % create the fundamental gaussian

n_groups = 40;
pq_map = zeros(3, n_groups^2);
temp = repmat(0:(n_groups-1), n_groups,1);          
pq_map(2,:) = temp(:);                              % assign p = 1,1,1,...2,2,2,...
pq_map(3,:) = repmat(0:(n_groups-1), 1,n_groups);   % assign q = 1,2,3...1,2,3...
[pq_map(1, :), I] = sort( pq_map(2,:) + pq_map(3,:) + 1, 'ascend' ); % group by asending p + q + 1
pq_map(2, :) = pq_map(2,I);
pq_map(3, :) = pq_map(3,I);
betas = n_0*k*sqrt( 1 - 2*Dn*(pq_map(1,:)/M) ); % Eq. 9

NMode = sum(betas > beta_low);
betas = betas(1:NMode);
pq_map = pq_map(:, 1:NMode);

%% calculation of mode profiles
E = zeros(N, N, NMode); % scalar field, independent of polarization
for ii = 1:NMode
    p = pq_map(2, ii);
    q = pq_map(3, ii);
    
    A = hermiteH(p,sqrt(2)*x/wo) .* hermiteH(q,sqrt(2)*y/wo) .* E0; % Eq. 4
    E(:,:,ii) = E_normalize(A);
end

%% construction of transmission matrix in LG representation with propagation constants
n_segments = numel(Length); % # of MMF segments
T = eye(2*n_count_modes); % consider both polarization

for pp = 1:n_segments                                                       % iterate over each piece of MMF segment
dl = Length(pp);
bend_r = Rho(pp);
theta = Theta(pp);
    if bend_r == inf                                                        % if there is no bending
        T = diag( exp(-1i * repmat(betas,1,2) * dl) ) * T;                  % use regular prop constants
    else
        C = GRIN_MMF_HG_coupling( k, n_0, wo, E, pq_map, bend_r );
        U = expm( (-1i*diag(repmat(betas,1,2)) + 1i*blkdiag(C,C))*dl );     % Eq. 24
        Rot = [cos(theta)*eye(n_count_modes), sin(theta)*eye(n_count_modes); -sin(theta)*eye(n_count_modes), cos(theta)*eye(n_count_modes)]; % Eq. 27
        M_proj = sectional_modal_proj(n_count_modes, pq_map, theta);        % Eq. 32
        T = (M_proj*Rot*U) * T;
    end
end

end



%% some self-defined functions
function E = E_normalize(E)
    E = E/sqrt(sum(abs(E(:)).^2));
end

function C = GRIN_MMF_HG_coupling(k, n_0, wo, E, pq_map, bend_r)
% this function implement Eq. 22 of
%   ??Shemirani et al., "Principal Modes in Graded-Index Multimode Fiber in Presence of Spatial- and Polarization-Mode Coupling"
%
% E is a N by N by n_count_modes array
% betas is a 1 by n_count_modes array
% dl is the segment length (unit: m)
% bend_r is the bending radius (unit: m)
% theta is the angle between x axis and the orientation of bending projected on x-y plane

n_count_modes = size(E, 3);
C = zeros(n_count_modes, n_count_modes);
coeff = 1i*k*n_0*(1/bend_r)*wo/2;

for ii = 1:n_count_modes % output
    for jj = 1:n_count_modes % input
        p = pq_map(2,ii);
        q = pq_map(3,ii);
        p_prime = pq_map(2,jj);
        q_prime = pq_map(3,jj);
        C(ii,jj) = coeff*( sqrt(p)*dirac(p,p_prime+1)*dirac(q,q_prime) + sqrt(p_prime)*dirac(p,p_prime-1)*dirac(q,q_prime) ); % coupling from jth mode to ith mode
    end
end
  
end

function c = dirac(a, b)
    if a == b
        c = 1;
    else
        c = 0;
    end
end

function K = sectional_modal_proj(n_count_modes, pq_map, theta)
% this function implement Eq. 30 of
%   ??Shemirani et al., "Principal Modes in Graded-Index Multimode Fiber in Presence of Spatial- and Polarization-Mode Coupling"
K = zeros(n_count_modes);
si = sin(theta);
co = cos(theta);
for ii = 1:n_count_modes % output
    p = pq_map(2,ii);
    q = pq_map(3,ii);
    for jj = 1:n_count_modes % input
        m = pq_map(2,jj);
        n = pq_map(3,jj);
        for k = 0:p
            for l = 0:q
                s = (k+q-l+m)/2;
                t = (p-k+l+n)/2;
                if (mod(s,1) == 0) && (mod(t,1) == 0) && s>=k && s>=(q-l) && s>=m && t>=p-k && t>=l && t>=n
                    K(ii,jj) = K(ii,jj) +...
                        (sqrt(factorial(p)*factorial(q)*factorial(m)*factorial(n))*((-1)^(p-k))*(co^k)*(si^(p-k))*(co^l)*(si^(q-l))) / ...
                        (factorial(s-k)*factorial(s-q+l)*factorial(s-m)*factorial(t-p+k)*factorial(t-l)*factorial(t-n));
                end
            end
        end
    end
end
    K = blkdiag(K,K);
end
