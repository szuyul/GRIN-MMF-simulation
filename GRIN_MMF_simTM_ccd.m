function T_record_basis = GRIN_MMF_simTM_ccd( lambda, D, Length, Rho, Theta, N, n_count_modes )
%GRIN_MMF_simTM_ccd simulates the complex images of a GRIN-MMF output pattern on a camera, subject to 
% the individual input realizations at the fiber input, which is a scanning focal grid. Only single
% polarization state is considered
%
% output_img = GRIN_MMF_simTM_HG( lambda, D, Length, Rho, Theta, N, n_count_modes )
% 
% output: 
% T_record_basis is a N^2 by N^2 2D transmission matrix in the recording basis
%
% input: 
% lambda is the wavelength (unit: m)
% D is the fiber core diameter of the MMF (unit: m)
% Length is an array of the segmental length of the MMF (unit: m)
% Rho is an array of the bending radius (unit: m) of each defined segment. If the segment is straight, Rho = inf
% Theta is an array of the orientation of the bending of the MMF (unit: rad.)
% N is the image dimension, e.g., 32
% n_count_modes is the number of counted modes for a GRIN MMF
% input_img is a N by N user defined 2D image
%
%
% 2020 Szu-Yu Lee
% BLCTO at Nokia Bell Labs

%% use eitehr Hermite Gaussian or Laugerre Gaussian basis
[ T, ~, ~, E, ~, ~ ]...
              = GRIN_MMF_simTM_HG( lambda, D, Length, Rho, Theta, N, n_count_modes );
%b = 2e-4;
%[ T, ~, ~, ~, E, ~, ~ ]...
%              = GRIN_MMF_simTM_LG( lambda, b, D, Length, Rho, N, n_count_modes );
   
%% transforming T into recording basis
T = T(1:n_count_modes, 1:n_count_modes);
E = reshape(E, [N^2, n_count_modes]);
T_record_basis = zeros(N^2, N^2);
for ii = 1:N^2
    temp = zeros(N,N);
    temp(ii) = 1;
    T_record_basis(:,ii) = E*T*E'*temp(:);
end

moviefixsc( reshape(T_record_basis, [N, N, N^2]) )

end

