function f = Calc_DrivingFocre(Phi,M,EPS,dx,dy,MethodFlag)

dFdPhi = Calc_dFdPhi(Phi);
Del2Phi = Calc_Del2Phi(Phi, dx,dy);

Mu = dFdPhi - 2 * EPS .* Del2Phi;

GradMu = Calc_Grad(Mu,dx,dy,MethodFlag);
GradM  = Calc_Grad(M,dx,dy,MethodFlag);
Div_GradMu = Calc_Div(GradMu,dx,dy,MethodFlag);
% g = M .* Div_GradMu + dot(GradM, GradMu,3) - 2*Reaction_rate;
f = M .* Div_GradMu + dot(GradM, GradMu,3);
