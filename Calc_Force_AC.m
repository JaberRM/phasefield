function f = Calc_Force_AC(Phi,C,M,EPS,dx,dy)


% G1 = (10*(C(:,:,1)-0.0)).^2 + 10;
% G2 = (10*(C(:,:,1)-1.0)).^2 + 10;
% % G = G1 * f(Phi) + G2 * (1-f(Phi))
% 
% dGdPhi = G1 .*  Calc_dfdPhi(Phi) -  G2 .* Calc_dfdPhi(Phi);

% DW
% G =  5*(phi-0.5).^4 - (phi-0.5).^2 ;
dGdPhi =  20 * (Phi-0.5).^3 - 2 * (Phi-0.5);


Del2Phi = Calc_Del2(Phi, dx,dy);

f = -M * (dGdPhi - 2 * EPS .* Del2Phi);
