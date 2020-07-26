function [ T, n_modes, lm_map, betas, E, img_size, refractive_idx ]...
              = GRIN_MMF_simTM_LG( lambda, b, D, Length, Rho, N, n_count_modes )
%GRIN_MMF_simTM_LG simulates the scalar LG modes, a set of Laguerre Gaussian functions, 
% of a GRIN MMF given the specifications and geometry. The LG solutions approximate well to propagation-invariant modes (PIMs.)
% The code is the implementation of theoretical GRIN MMF model, which can be found in:
% Flaes et al., "Robustness of Light-Transport Processes to Bending Deformations in Graded-Index Multimode Waveguides, "
% Phys. Rev. Lett. 120, 233901 (2018).
% The fiber transmission of a deformed MMF is simulated based on a perturbation theory by correcting the propagation constants, 
% while the LG mode profiles do not change much.
%
% [ T, NMode, lm_map, betas, E, img_size, refractive_idx ]...
%            = GRIN_MMF_simTM_LG( lambda, D, NA, Length, Rho, Theta, N, n_count_modes ) 
% 
% outputs:
% T is the transmission matrix of a straight or deformed GRIN MMF in LG representation at a certain wavelength
% n_modes is the number of total LG modes
% lm_map is a 3 by n_modes matrix, with each column as the characteristics of one particular mode
%                                      1st row as the angular momentum, 2m + |l|
%                                      2nd row as the orbital angular momentum
%                                      3rd row as the # of radial nodes 
% beta is a 1 by n_modes matrix, the propagation constants (unit: m) of each LG mode
% E is a N by N by n_modes matrix, the transverse electric scalar field distribution of each LG mode in one polarization
% img_size is the physical size (unit: m) of the transverse field images
% refractive_idx is the refractive index profile over the img_size
%
% inputs:
% lambda is the operating wavelength (unit: m)
% b is the parabolic parameter of the MMF refractive index 
% D is the fiber core diameter of the MMF (unit: m)
% Length is an array of the segmental length of the MMF (unit: m)
% Rho is an array of the bending radius (unit: m) of each defined segment. If the segment is straight, Rho = inf
% Theta is an array of the orientation of the bending of the MMF (unit: rad.)
% N is the image dimension
% n_count_modes is the number of counted modes for a GRIN MMF
%
%
% 2020 Szu-Yu Lee
% BLCTO at Nokia Bell Labs

%% GRIN MMF parameter initialization (unit in m)
k = 2*pi/lambda;
n_0 = 1.48; % central refractive index 

alpha = k*n_0/b;
sig = 1; % R(+1) or L(-1) circular pol.

img_size = 1*D;
range = linspace(-img_size/2, img_size/2, N);
[x, y] = meshgrid(range);
refractive_idx = n_0*(1 - range.^2/b^2);
[theta, rho] = cart2pol(x,y);

lmax = 20;
mmax = 15;
[l_grid, m_grid] = meshgrid(-lmax:lmax, 0:mmax);
lm_map = zeros(3, numel(l_grid));
[lm_map(1, :), I] = sort( abs(l_grid(:)) + 2*m_grid(:), 'ascend' ); % group by asending |l| + 2*m
lm_map(2, :) = l_grid(I);
lm_map(3, :) = m_grid(I);
lm_map = lm_map(:, 1:n_count_modes);

%% calculate the scalar field of each mode
n_modes = size(lm_map, 2);
E = zeros(N, N, n_modes);

for ii = 1:n_modes
    l = lm_map(2, ii); 
    m = lm_map(3, ii); 
    
    E(:,:,ii) = GRIN_MMF_LG_fields( alpha, l, m, rho, theta );
end
betas = GRIN_MMF_LG_betas( lambda, n_0, b, alpha, lm_map, sig, inf ); % propagation constants for straight MMF

%% construction of transmission matrix in LG representation with propagation constants
%mode_to_x = reshape(E, [N*N, n_modes]);

n_segments = numel(Length); % # of MMF segments
T = eye(n_modes);

for pp = 1:n_segments                                                       % iterate over each piece of MMF segment
    if Rho(pp) == inf                                                       % if there is no bending
        T = diag( exp(-1i * betas * Length(pp)) ) * T;                       % use regular prop constants
    else
        bend_r = Rho(pp);
        [~, betas_prime] = GRIN_MMF_LG_betas( lambda, n_0, b, alpha, lm_map, sig, bend_r );
        T = diag( exp(-1i * betas_prime * Length(pp)) ) * T;
    end
end
    
end



%% some self-defined functions
function E = GRIN_MMF_LG_fields( alpha, l, m, r, theta )
    ar_2 = alpha*(r.^2);

    c1 = sqrt( (alpha/(2*pi))*( (2*factorial(m))/factorial(m+abs(l)) ) );
    c2 = exp( -ar_2/2 );
    c3 = ar_2.^(abs(l)/2);
    c4 = laguerreL(m, abs(l), ar_2); % Generalized Laguerre polynomials (https://en.wikipedia.org/wiki/Laguerre_polynomials)
    c5 = exp(1i*l*theta);
    E = c1*c2.*c3.*c4.*c5; % scalar field
    E = E_normalize(E);
end

function [betas, betas_prime] = GRIN_MMF_LG_betas( lambda, n_0, b, alpha, lm_map, sig, bend_r )
% lm_map is a 3 by n_modes matrix, with each column as the characteristics of one particular mode
%                                      1st row as the angular momentum, 2m + |l|
%                                      2nd row as the orbital angular momentum
%                                      3rd row as the # of radial nodes
    k = 2*pi/lambda;
    l_s = lm_map(2,:);
    m_s = lm_map(3,:);
    
    betas = sqrt((k*n_0)^2 - 2*alpha*(abs(l_s) + 2*m_s + 1)); % 1 by n_modes matrix, propagation constants in straight MMF
    curv = 1/bend_r; % curvature
    SO_corrections = -(l_s*sig + 1)/(2*k*n_0*b^2); % 1 by n_modes matrix, spin-orbital interactions 
    betas_prime = betas + (curv^2)*( (k*n_0*b^2)/2 - 9*b*(abs(l_s) + 2*m_s + 1)/4 ) + SO_corrections; % 1 by n_modes matrix, propagation constants in bent MMF
end

function E = E_normalize(E)
    E = E/sqrt(sum(abs(E(:)).^2));
end