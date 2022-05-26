clc;  clearvars;
close all;

restoredefaultpath

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% define the geometry of sources and the target region

%%% filamentary sources: a total of N_source sources placed on a circle of
%%% radius rho_source
R0 = 1;
Z0 = 0;
N_source = 11;

rho_source = .6;

theta_source = linspace(0,2*pi,N_source+1).';
theta_source = theta_source(1:end-1);

RZ_source = [R0 Z0] + rho_source*[cos(theta_source) sin(theta_source)];

% chose sinusoidal test current
I_source = cos(theta_source);


%%% target region: circle of radius rho_target
rho_target = .9*rho_source;

theta_target = linspace(0,2*pi,100).';
% % theta_target = theta_target(1:end-1);

RZ_target = [R0 Z0] + rho_target*[cos(theta_target) sin(theta_target)];


%%% Plot geometry
ind_fig_1 = floor(1e+5*rand);
figure(ind_fig_1)
plot(RZ_target(:,1),RZ_target(:,2),'.-r')
hold on; axis equal
plot(RZ_source(:,1),RZ_source(:,2).','ok', 'LineWidth',2)
xlabel('$r$ [m]')
xlabel('$z$ [m]')
legend('target region', 'sources')



%% Compute magnetic vector potential (phi component) inside the target region
% define grid 
[RR,ZZ] = meshgrid(linspace(min(RZ_target(:,1)),max(RZ_target(:,1)),50), ...
    linspace(min(RZ_target(:,2)),max(RZ_target(:,2)),50));

vec_RR = RR(:);
vec_ZZ = ZZ(:);


%%% use routine in /fun_Green_Fortran_ver_1.01
addpath('./fun_Green_Fortran_ver_1.01')

OPT_parallel = 1; % 1->parallel, 0->serial
n_threads = 8;    % no. of threads (effectively considered only of OPT_parallel=1)

Aphi_inside = fun_Green_filament_Aphi_SP_f90(length(RZ_source(:,1)),...  % no. of source points
    RZ_source(:,1),... % r of source points
    RZ_source(:,2),... % z of source points
    I_source,... % currents of source points
    length(vec_RR),... % no. of target points
    vec_RR,... % no. of source points
    vec_ZZ,... % no. of source points
    OPT_parallel,... 
    n_threads); 

% replace the solution with NaN for grid points outside the target domain
Aphi_inside(~inpolygon(vec_RR,vec_ZZ,RZ_target(:,1),RZ_target(:,2))) = NaN;


%%% Plot result
nsurf_cont = 30;
figure(ind_fig_1)
contourf(RR,ZZ,reshape(Aphi_inside,size(RR)),nsurf_cont)
colormap(jet(nsurf_cont))
colorbar vert



%% Compute magnetic vector potential on the boundary of the target region

Aphi_boundary = fun_Green_filament_Aphi_SP_f90(length(RZ_source(:,1)),...  % no. of source points
    RZ_source(:,1),... % r of source points
    RZ_source(:,2),... % z of source points
    I_source,... % currents of source points
    size(RZ_target,1),... % no. of target points
    RZ_target(:,1),... % no. of source points
    RZ_target(:,2),... % no. of source points
    OPT_parallel,... 
    n_threads); 


%%%
ind_fig_2 = floor(1e+5*rand);
figure(ind_fig_2)
plot(theta_target,Aphi_boundary, 'k', 'LineWidth', 2)
hold on
xlabel('$\theta$')
ylabel('$A_{\phi}$')


%% DFT based filter along theta anfgle

%%% use routine in /routines_Fourier_v_1.0
addpath('./routines_Fourier_v_1.0')

mm = [-20:20]; %harmonics
[c_n]=fun_FFT(Aphi_boundary,mm.',theta_target); % Fourier harmonics
Aphi_boundary_filt = fun_FFT_inv(c_n,mm,theta_target); % filtered physical quantity

figure(ind_fig_2)
plot(theta_target,Aphi_boundary_filt, '--r', 'LineWidth', 2)
legend('vector pot. from Biot-Savart', sprintf('vector pot. filtered with %i harmonics',length(mm)))


%% Create a large dataset for training the Feed-Forward Neural Net (FFNN)


N_source_all = 1:10;
rho_source_all = linspace(1.1*rho_target,2*rho_target,10);
freq_all = 1:10;
mm = [-20:20]; %harmonics

dadaset_A_boundary = zeros(length(mm),length(N_source_all)*length(rho_source_all)*length(freq_all));
dadaset_A_inside = zeros(length(vec_RR),length(N_source_all)*length(rho_source_all)*length(freq_all));

tic
fprintf('Assembling the database ... \n')

kk = 1;
for ii = 1:length(N_source_all)
    fprintf('done %i of %i \n',ii*length(rho_source_all)*length(freq_all), ...
        length(N_source_all)*length(rho_source_all)*length(freq_all))
    for jj = 1:length(rho_source_all)
        for hh = 1:length(freq_all)

            % explore all combinations
            theta_source = linspace(0,2*pi,N_source_all(ii)+1).';
            theta_source = theta_source(1:end-1);
            RZ_source = [R0 Z0] + rho_source_all(jj)*[cos(theta_source) sin(theta_source)];
            I_source = cos(freq_all(hh)*theta_source);

            %%% Aphi inside the target domain
            Aphi_inside = fun_Green_filament_Aphi_SP_f90(length(RZ_source(:,1)),...  % no. of source points
                RZ_source(:,1),... % r of source points
                RZ_source(:,2),... % z of source points
                I_source,... % currents of source points
                length(vec_RR),... % no. of target points
                vec_RR,... % no. of source points
                vec_ZZ,... % no. of source points
                1,...
                n_threads);

            % replace the solution with NaN for grid points outside the target domain
            Aphi_inside(~inpolygon(vec_RR,vec_ZZ,RZ_target(:,1),RZ_target(:,2))) = NaN;


            %%% Aphi on boundary
            Aphi_boundary = fun_Green_filament_Aphi_SP_f90(length(RZ_source(:,1)),...  % no. of source points
                RZ_source(:,1),... % r of source points
                RZ_source(:,2),... % z of source points
                I_source,... % currents of source points
                size(RZ_target,1),... % no. of target points
                RZ_target(:,1),... % no. of source points
                RZ_target(:,2),... % no. of source points
                0,...
                n_threads);

            if 0
                [c_n]=fun_FFT(Aphi_boundary,mm.',theta_target); % Fourier harmonics
                Aphi_boundary_filt = fun_FFT_inv(c_n,mm,theta_target); % filtered physical quantity
                figure
                plot(theta_target,Aphi_boundary, 'k', 'LineWidth', 2)
                hold on
                xlabel('$\theta$')
                ylabel('$A_{\phi}$')
                plot(theta_target,Aphi_boundary_filt, '--r', 'LineWidth', 2)
                legend('vector pot. from Biot-Savart', sprintf('vector pot. filtered with %i harmonics',length(mm)))
            end

            %%%
            dadaset_A_boundary(:,kk) = c_n;
            dadaset_A_inside(:,kk) = Aphi_inside;

            kk = kk+1;

        end
    end
end
toc








