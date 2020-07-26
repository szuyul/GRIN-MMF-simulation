%% simulating LG modes in GRIN MMF
lambda = 1.54e-6;
b = 4e-4; % parabolic parameter
D = 8e-5; % diameter (https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=358)

N = 30;
Length = 1;
Rho = inf;
% for a bent MMF
%Length = 0.1*ones(1,10);
%Rho = 0.2 + 0.05*rand(1,numel(Length));

n_count_modes = 45;
[ T, n_modes, lm_map, propconst, E, img_size, refractive_idx ]...
              = GRIN_MMF_simTM_LG( lambda, b, D, Length, Rho, N, n_count_modes );

%% showing the refractive index
range = linspace(-img_size/2, img_size/2, N);
lmax = max(lm_map(2,:));
qmax = max(lm_map(3,:));

close all
figure
plot(range*1e6, refractive_idx)
title('refractive index profile')
xlabel('r (\mum)')
grid on

%% show all the sorted LG modes at the same time
figure('Position', [100, 500, 1100, 400]);
for ii = 1:n_modes
    temp = E(:,:,ii);                                                       % choose a scalar field to plot
    
    figure(2)
    pind = (lm_map(3,ii))*(2*lmax+1) + lm_map(2,ii) + lmax + 1;
    subplot(qmax+1, 2*lmax+1, pind);
    complex_imagesc( temp );
    title(['LG l = ',num2str(lm_map(2,ii)),', m = ',num2str(lm_map(3,ii))]);
    axis off
end

%% show the sorted LG modes sequentially 
figure
moviefixsc( E )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% simulating HG modes in GRIN MMF
lambda = 1.55e-6;
D = 50e-6; 
Length = 0.02*ones(1, 50);
Rho = 0.8 + .1*(rand(1, numel(Length))-0.5);
Theta = 2*pi*(rand(1,numel(Length))-0.5);
N = 30;
n_count_modes = 55;

[ T, pq_map, betas, E, img_size, refractive_idx ]...
              = GRIN_MMF_simTM_HG( lambda, D, Length, Rho, Theta, N, n_count_modes );
pmax = max(pq_map(2,:));
qmax = max(pq_map(3,:));

%% show all the sorted HG modes at the same time
figure('Position', [100, 500, 1100, 400]);
for ii = 1:n_count_modes
    temp = E(:,:,ii);                                                       % choose a scalar field to plot
    
    pind = pq_map(2,ii)*(qmax+1) + pq_map(3,ii) + 1;
    subplot(pmax+1, qmax+1, pind);
    complex_imagesc( temp );
    title(['HG p = ',num2str(pq_map(2,ii)),', q = ',num2str(pq_map(3,ii))]);
    axis off
end

%% show the sorted HG modes sequentially 
figure
moviefixsc( E )
