clc; close all; clearvars;
system('rm -rf ./Pics');
system('mkdir ./Pics');

%% System configuration
Scale = 0.2;
Im_ice = imread('map_ice.png');
Im_ice = im2bw(Im_ice,0.5);
Im_ice = imresize(Im_ice,Scale);
% imshow(Im_ice)

Im_land = imread('map_land.png');
Im_land = im2bw(Im_land,0.5);
Im_land = imresize(Im_land,Scale);
% imshow(Im_land)

Land = Im_land == 0 ;
Ice = Im_ice == 1;

[Lx,Ly] = size(Im_land);
Lx = Lx - 1;
Ly = Ly - 1;

Param.dx = 1;
Param.dy = 1;

% Simulation time
Time = 4000000000;      % s
Param.dt = 10000 ;  % s

% Number of conservative (Cahn-Hilliard) fields
NCH = 0;

% Number of non-conservative (Allen-Cahn) fields
NAC = 1;

% Number of thermal fields
NTH = 1;

% Gradient coeff.
EPS_CH(1:NCH) = 2e0;
EPS_AC(1:NAC) = 1e1;

% Mobility
M0_CH(1:NCH) = 5e-1;
M0_AC(1:NAC) = 1e-7;

% Output options
EveryNStep = 200;

% Thermal
TH.Cp_water = 4.187e3; % heat capacity of water J/kgK;
TH.Cp_ice = 2.108e3;   % heat capacity of ice J/kgK
TH.Cp_land = 8.0e2;    % heat capacity of land J/kgK

TH.Q =  334000;        % latent heat of ice J/Kg

TH.K_water = 0.5918;    % thermal conductivity of water W·m−1·K−1
TH.K_ice =   2.0914;    % thermal conductivity of water W·m−1·K−1
TH.K_land =  0.255;     % thermal conductivity of water W·m−1·K−1

TH.Rho_water = 997;  %kg/m3
TH.Rho_ice =  916.8; %Kg/m3
TH.Rho_land = 2000; %Kg/m3

%% Initialize variables
X0 = 0 : Param.dx : Lx;
Y0 = 0 : Param.dx : Ly;

Nx = length(X0);
Ny = length(Y0);

X = meshgrid(X0,Y0);
Y = meshgrid(Y0,X0)';

% C = zeros(Nx, Ny , NCH)  + 0.5 + (rand(Nx,Ny,NCH)-0.5)/100;

Phi = zeros(Nx, Ny , NAC);
T = ones(Nx, Ny , NTH) * 273;

Flag_Ice = repmat(Ice,1,1,NTH);
Flag_Land= repmat(Land,1,1,NTH);

[A,B] = gradient(Land);
Edge = A.^2 + B.^2;

Phi(Flag_Ice) = 1;
T(Flag_Ice) = 273-40;

%% Run
Iter_AC = 0;
Iter_CH = 0;
Iter_TH = 0;

Tol_AC = 0;
Tol_CH = 0;
Tol_TH = 0;

Step = 0;
cnt2 = 0;
CT=cbrewer('div','PiYG',500);
set(gcf,'position',[ 214   105   861   693]);


for time = 0: Param.dt : Time
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
        [Phi(:,:,i), Tol_AC,Iter_AC, converged] = Update_AC( Phi(:,:,i), T , M0_AC(i), EPS_AC(i), Param );
        if converged == 0
            fprintf('Allen-Cahn diverged! Iter_AC=%i , Tol_AC=%g \n', Iter_AC , Tol_AC);
            break;
        end
    end
    
    for i = 1:NTH
        [T(:,:,i), Tol_TH,Iter_TH, converged] = Update_TH( T(:,:,i), Phi, Param, TH, Flag_Land);
        if converged == 0
            fprintf('Allen-Cahn diverged! Iter_TH=%i , Tol_TH=%g \n', Iter_TH , Tol_TH);
            break;
        end
    end
    
    if converged == 0
        fprintf('diverged! \n');
        break;
    end
    
    % postprocessing & output
    if mod(Step,EveryNStep)==1
        cnt2 = cnt2 + 1;
        fprintf('Step = %i , time = %2.3e ,  Iter_AC = %i , Iter_TH = %i , Tol_AC = %2.3e , Tol_TH = %2.3e \n', Step, time, Iter_AC , Iter_TH , Tol_AC , Tol_TH );
        cla;

        subplot(1,2,1)
        hold on
        
%         surf(X,Y,1000*Edge','edgecolor','none'); view(2); xlim([0,max(X(:))]); ylim([0,max(Y(:))]); daspect([1 1 100]); 
        Id = Edge' > 0.5;
        scatter3(X(Id),Y(Id),Edge(Id)+1000,'.','Linewidth',0.1,'MarkerEdgeColor',[0 0 0])
        
        surf(X,Y,T(:,:,1)','edgecolor','none'); view(2); xlim([0,max(X(:))]); ylim([0,max(Y(:))]); daspect([1 1 100]);
        colormap(gca,'parula');
%         caxis([ min(T(:)), max(T(:))])
        xticks([]); xlabel('');
        yticks([]); ylabel('');
        
        subplot(1,2,2)
        hold on
        surf(X,Y,Phi(:,:,1)','edgecolor','none'); view(2); xlim([0,max(X(:))]); ylim([0,max(Y(:))]); daspect([1 1 100]);
        scatter3(X(Id),Y(Id),Edge(Id)+1000,'.' ,'Linewidth',0.1,'MarkerEdgeColor',[0 0 0])
        colormap(gca,'parula');
        xticks([]); xlabel('');
        yticks([]); ylabel('');
        
        drawnow



    end
    
end


