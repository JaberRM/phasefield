function f = Calc_Force_TH(T,K,dx,dy)

GradT = Calc_Grad(T,dx,dy,1);
GradK  = Calc_Grad(K,dx,dy,1);
Div_GradT = Calc_Div(GradT,dx,dy,-1);
% f = M .* Div_GradMu + dot(GradK, GradT,3) - 2*Reaction_rate;
% f = K.* Div_GradT + dot(GradK, GradT,3);
f = K .* Div_GradT ;
