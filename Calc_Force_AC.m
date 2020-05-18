function f = Calc_Force_AC(Phi,x,theta,Param)


GradPhi = Calc_Grad(Phi,Param.dx,Param.dy,1);
dPhidx = GradPhi(:,:,1);
dPhidy = GradPhi(:,:,2);

Temp = (Param.XOp - Param.XOmin + Param.eta) ./ (x + Param.XOmin + Param.eta); % Eq 10
nu0 = 4*Param.Xi * ( Temp.^12 -  Temp.^6 ); % Eq 10
       
m = Param.alpha/pi  .* Phi .* atan(nu0) + 0.5 * (Phi-1); % Eq 9

sigma = 1 + Param.delta * cos( Param.J * (theta-Param.theta0) ); % Eq 6 
epsilon = Param.eps_bar * sigma; % Eq 5
epsilon_p = - Param.eps_bar * Param.delta * sin( Param.J * (theta-Param.theta0) ) * Param.J;


% force based on Eq 4 (Terms are the right hand side)
Term1 = Calc_Grad(epsilon .* epsilon_p .* dPhidy, Param.dx,Param.dy,1);
Term1 = -Term1(:,:,1);

Term2 = Calc_Grad(epsilon .* epsilon_p .* dPhidx,Param.dx,Param.dy,1);
Term2 = Term2(:,:,2);

Temp = epsilon.^2 .* GradPhi;
Term3 = Calc_Div(Temp,Param.dx,Param.dy,1);

Term4 = Phi .* (1-Phi) .* (Phi - 1/2 + m);

% force
f = 1/Param.tau * (Term1 + Term2 + Term3 + Term4);
