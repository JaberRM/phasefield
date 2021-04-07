clc; close all; clearvars;
system('rm -rf ./Pics3');
system('mkdir ./Pics3');

%% System configuration
% Simulation Geometry
Lx = 200;
Ly = 200;

Param.dx = 1;
Param.dy = 1;

% Simulation time
T = 4000;      % s
Param.dt = 1e-2 ;  % s

% Number of conservative (Cahn-Hilliard) fields
NCH = 0;

% Number of non-conservative (Allen-Cahn) fields
NAC = 1;

% Gradient coeff.
EPS_CH(1:NCH) = 2e0;
EPS_AC(1:NAC) = 1e0;

% Mobility
M0_CH(1:NCH) = 5e-1;
M0_AC(1:NAC) = 5e-1;

% Output options
EveryNStep = 1000;

% Mech
% Young = 1000;
% nu = 0.2;
% Param.C = [1-nu     nu      0;
%     nu       1-nu    0;
%     0        0       1-2*nu] * Young / ((1+nu) * (1-2*nu));

%% Initialize variables
X0 = 0 : Param.dx : Lx;
Y0 = 0 : Param.dx : Ly;

Nx = length(X0);
Ny = length(Y0);

X = meshgrid(X0,Y0);
Y = meshgrid(Y0,X0)';

C = zeros(Nx, Ny , NCH)  + 0.5 + (rand(Nx,Ny,NCH)-0.5)/100;
Phi = zeros(Nx, Ny , NAC) + 0.5 + (rand(Nx,Ny,NAC)-0.5)/100;


%% Run
Iter_AC = 0;
Iter_CH = 0;
Tol_AC = 0;
Tol_CH = 0;
Step = 0;
cnt2 = 0;
CT=cbrewer('div','PiYG',500);
set(gcf,'position',[ 214   105   861   693]);


for time = 0: Param.dt : T
    Step = Step + 1;
    
    for i = 1:NCH
        [C(:,:,i), Tol_CH, Iter_CH, converged]= Update_CH( C(:,:,i) , Phi , M0_CH(i), EPS_CH(i), Param );
        if converged == 0
            fprintf('Cahn-Hilliard diverged! Iter_CH=%i , Tol_CH=%g \n', Iter_CH , Tol_CH);
            break;
        end
    end
    %
    for i = 1:NAC
        [Phi(:,:,i), Tol_AC,Iter_AC, converged] = Update_AC( Phi(:,:,i), C , M0_AC(i), EPS_AC(i), Param );
        if converged == 0
            fprintf('Allen-Cahn diverged! Iter_AC=%i , Tol_AC=%g \n', Iter_AC , Tol_AC);
            break;
        end
    end
    
    if converged == 0
        fprintf('Mech diverged! Iter_Mech=%i , Tol_Mech=%g \n', Iter_Mech , Tol_Mech);
        break;
    end
    
    % postprocessing & output
    if mod(Step,EveryNStep)==1
        cnt2 = cnt2 + 1;
        fprintf('Step = %i , time = %2.3e ,  Iter_AC = %i , Iter_CH = %i , Tol_AC = %2.3e , Tol_CH = %2.3e \n', Step, time, Iter_AC , Iter_CH , Tol_AC , Tol_CH );
        %         fprintf('Step = %i , time = %2.3e ,  Iter_Mech = %i , Tol_Mech = %2.3e \n', Step, time, Iter_Mech , Tol_Mech );
        cla;
        surf(X,Y,Phi(:,:,1)','edgecolor','none'); view(2); xlim([0,max(X(:))]); ylim([0,max(Y(:))]); daspect([1 1 100]);
        %         colorbar;
        %         colormap(gca,'parula');
        colormap(gca, CT);
        xticks([]); xlabel('');
        yticks([]); ylabel('');
        drawnow
        saveas(gcf,sprintf('./Pics3/%i.jpg',cnt2))
    end
    
end


