clc; close all; clearvars;

%% System configuration
% Simulation Geometry
Lx = 50e-6;
Ly = 50e-6;

dx = 1e-6;
dy = 1e-6;

% Simulation time
T = 10;      % s
dt = 1e-4 ;  % s

% Number of Cahn-Hilliard fields
NCH = 1;

% Gradient coeff.
EPS(1:NCH) = 1e-20;

% Mobility
M0(1:NCH) = 1e-6;

% Output options
EveryNStep = 1000;

%% Initialize variables
X = 0 : dx : Lx;
Y = 0 : dx : Ly;

Nx = length(X);
Ny = length(Y);

Phi = zeros(Nx, Ny , NCH) + rand(Nx,Ny,NCH)/100;
M = ones(Nx,Ny,NCH)*1e-4;

%% Run
Step = 0;
for time = 0: dt : T
    Step = Step + 1;

    for i = 1:NCH
        Phi(:,:,i) = Update_CH( Phi(:,:,i) , M(:,:,i), EPS(i), dx , dy , dt , 0);
    end
    
    % postprocessing & output
    if mod(Step,EveryNStep)==1
        fprintf('Step = %i,  time = %g\n', Step, time);
        cla;
        clf;
        surf(X,Y,Phi(:,:,1),'edgecolor','none'); view(2); xlim([0,max(X)]); ylim([0,max(Y)]); colorbar;  axis equal; colormap(gca,'parula');  drawnow
    end
    
end


