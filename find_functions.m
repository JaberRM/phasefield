clc
close all
clearvars


%% Gibbs functions for each phase
C = 0 : 0.001:1;
G1 = (5*(C-0.0)).^2 + 10;
G2 = (5*(C-1)).^2 + 10;

figure; hold on
plot(C,G1);
plot(C,G2);

%% indicator function for phase transformation
% Gibbs = G_1(c,T,..) * f(Phi) + G_2(c,T,..) * (1-f(Phi)) + 0.5 K_C grad(C)
% + 0.5 K_Phi grad(phi)^2 + h(phi)
% interface energy daniel schneider nestler ingo
% KKS
% f = phi^2 * ( a*phi^2 + b*phi + c)

a = 3;
c = a + 3;
b = 1 - a - c;

phi = -0.5:0.01:1.5;


f = phi.^2 .* (a*phi.^2 + b* phi + c);
% y2 = (1-phi).^2 .* (a*(1-phi).^2 + b* (1-phi) + c);
% y3 = G1(C==0.8) * y1 + G2(C==0.8)*y2;


% f = phi.^2 .* (3-2*phi);
% y2 = 1-y1;

c0 = 0.4;
g3 = G1(C==c0) * f + G2(C==c0)*(1-f);


figure
hold on
plot(phi,f);
% plot(phi,y2);
plot(phi,g3);

xlim([-1,2])











