clc; close all; clearvars;

%% System configuration
% Simulation Geometry
Lx = 50;
Ly = 50;

dx = 1;
dy = 1;

% Simulation time
T = 400;      % s
dt = 1e-2 ;  % s

% Number of conservative (Cahn-Hilliard) fields
NCH = 1;

% Number of non-conservative (Allen-Cahn) fields
NAC = 1;

% Gradient coeff.
EPS_CH(1:NCH) = 1e1;
EPS_AC(1:NAC) = 1e1;

% Mobility
M0_CH(1:NCH) = 1e-1;
M0_AC(1:NAC) = 1e-3;

% Output options
EveryNStep = 100;

%% Initialize variables
X = 0 : dx : Lx;
Y = 0 : dx : Ly;

Nx = length(X);
Ny = length(Y);

C = zeros(Nx, Ny , NCH) + rand(Nx,Ny,NCH)/100 + 0.45;
Phi = zeros(Nx, Ny , NCH) + rand(Nx,Ny,NCH)/100 + 0.01;

%% Run
Step = 0;
for time = 0: dt : T
    Step = Step + 1;
    
    for i = 1:NCH
        [C(:,:,i), Tol, Iter]= Update_CH( C(:,:,i) , Phi , M0_CH(i), EPS_CH(i), dx , dy , dt , 0);
    end
    
    for i = 1:NAC
        [Phi(:,:,i), Tol,Iter] = Update_AC( Phi(:,:,i), C , M0_AC(i), EPS_AC(i), dx , dy , dt );
    end
    
    % postprocessing & output
    if mod(Step,EveryNStep)==1
        fprintf('Step = %i , time = %g\n', Step, time);
        clf;
        subplot(1,2,1)
        cla;
        surf(X,Y,C(:,:,1),'edgecolor','none'); view(2); xlim([0,max(X)]); ylim([0,max(Y)]); daspect([1 1 100]); colorbar; colormap(gca,'parula');
        subplot(1,2,2)
        cla;
        surf(X,Y,Phi(:,:,1),'edgecolor','none'); view(2); xlim([0,max(X)]); ylim([0,max(Y)]); daspect([1 1 100]); colorbar; colormap(gca,'parula');
        drawnow
    end
    
end


