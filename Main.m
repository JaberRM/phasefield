% implementation of the model in 
% Lu, F., Wen, L., Li, J., Wei, J., Xu, J., & Zhang, S. (2016). 
% Numerical simulation of iron whisker growth with changing oxygen content 
% in iron oxide using phase-field method. 
% Computational Materials Science, 125, 263–270. 
% https://doi.org/10.1016/j.commatsci.2016.09.003

clc; close all; clearvars;
%% System configuration
% Parameters in table 2
Param.tau = 0.0003;
Param.theta0 = pi/4;
Param.J = 4;
Param.delta = 0.1;
Param.alpha = 0.9;
% Param.alpha = 1.9;
Param.Xi = 2.450;
Param.XOp = 0.023;
Param.eta = 0.710;
Param.D_s = 1e-9;
Param.D_int = 1e-5;

Param.D_s = 1e-4;
Param.D_int = 1e-0;

% Param.D_s = 1e-3;
% Param.D_int = 1e0;

% Other parameters from the paper
% Param.dt = 0.0002 ;  % s
Param.dt   = 0.000001 ;  % s

Lx = 3;
Ly = 3;

Param.dx = 0.01;
Param.dy = 0.01;

Param.eps_bar = 0.02;
% Param.eps_bar = 0.04;

% Param.K = 7.3e-3;
Param.K = 7.3e-1;

%unknown params
Param.XOmin  = 0.002;

% Simulation time
T = 2;
% Output options

EveryNStep = 50;

%% Initialize variables
X0 = -Lx/2 : Param.dx : Lx/2;
Y0 = -Ly/2 : Param.dx : Ly/2;

Nx = length(X0);
Ny = length(Y0);

X = meshgrid(X0,Y0);
Y = meshgrid(Y0,X0)';


InCirc = X.^2 + Y.^2 <= 1;

x   = zeros(Nx, Ny);
Phi = zeros(Nx, Ny);
Phi(InCirc) = 1;
x(InCirc) = 0.6;

%% Run
Step = 0;
for time = 0: Param.dt : T
    Step = Step + 1;
    
    Phi0 = Phi;
    
    [Phi, Tol_AC, Iter_AC, converged] = Update_AC( Phi , x , Param );
    if converged == 0
        fprintf('Allen-Cahn diverged! Iter_AC=%i , Tol_AC=%g \n', Iter_AC , Tol_AC);
        break;
    end
    
    dPhidt = (Phi-Phi0) / Param.dt;

    [x,   Tol_DF, Iter_DF, converged] = Update_DF( x , Phi, dPhidt , Param );
    if converged == 0
        fprintf('Diffusion diverged! Iter_DF=%i , Tol_DF=%g \n', Iter_DF, Tol_DF);
        break;
    end
    
    
    
    % postprocessing & output
    if mod(Step,EveryNStep)==0
        fprintf('Step = %i , time = %2.3e ,  Iter_AC = %i , Iter_DF = %i , Tol_AC = %2.3e , Tol_DF = %2.3e \n', Step, time, Iter_AC , Iter_DF , Tol_AC , Tol_DF );
        clf;
        subplot(3,1,1)
        cla;
        surf(X,Y,Phi,'edgecolor','none'); view(2); xlim([min(X0),max(X0)]); ylim([min(Y0),max(Y0)]); daspect([1 1 100]); colorbar; colormap(gca,'jet'); set(gca,'fontsize', 15); ylabel('\phi');
        subplot(3,1,2)
        cla;
        surf(X,Y,x,'edgecolor','none'); view(2); xlim([min(X0),max(X0)]); ylim([min(Y0),max(Y0)]); daspect([1 1 100]); colorbar; colormap(gca,'jet'); set(gca,'fontsize', 15); ylabel('x_o');
        
        subplot(3,1,3)
        cla;
        
        Temp = (Param.XOp - Param.XOmin + Param.eta) ./ (x + Param.XOmin + Param.eta); % Eq 10
        nu0 = 4*Param.Xi * ( (Temp).^12 - ( (Temp).^6) ); % Eq 10
        m = Phi .* Param.alpha/pi .* atan(nu0) + 0.5 * (Phi-1) +0.2 ; % Eq 9

        surf(X,Y,m,'edgecolor','none'); view(2); xlim([min(X0),max(X0)]); ylim([min(Y0),max(Y0)]); daspect([1 1 100]); colorbar; colormap(gca,'jet'); set(gca,'fontsize', 15); ylabel('m');
        set(gcf,'position',[   440     1   438   797]);
        drawnow
    end
    
end


