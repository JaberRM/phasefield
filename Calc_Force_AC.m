function f = Calc_Force_AC(Phi,T,M,EPS,dx,dy,  Mat, Flag_Land)

% https://www.researchgate.net/publication/228365777_Phase-Field_Simulation_of_Solidification_1
% G = WA * g(phi) + LA (TM - T) / TM * p(phi)
% g(phi) = phi^2 (1-Phi)^2
% p(phi) = phi^3 * (6*phi^2 - 15Phi + 10)
WA = 1e2;

dgdphi = 2*Phi .* (1-Phi).^2 - 2 * Phi.^2 .* (1-Phi);
dpdphi = 3*Phi.^2 .*  (6*Phi.^2 - 15*Phi + 10) +  Phi.^3 .* (12*Phi - 15);
dGdPhi = WA * dgdphi + Mat.Q * (Mat.Tm - T) / Mat.Tm .* dpdphi;

Del2Phi = Calc_Del2(Phi, dx,dy);

f = - M * (dGdPhi - 2 * EPS .* Del2Phi);
