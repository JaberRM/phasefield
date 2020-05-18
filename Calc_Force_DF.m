function f = Calc_Force_DF(x,Phi,Param)

% Diffusivity based on Lu et al.
D = zeros(size(x));
D(Phi>=0.95) = Param.D_s;
D(Phi<0.95 & Phi>0.05) = Param.D_int;
D(Phi<=0.05) = 0;

% f = Del2(Dx) from Eq 11 (this is wrong!)
% Del2Dx = Calc_Del2(D.*x, Param.dx,Param.dy);
% f =  Del2Dx ; %Eq 11 this is wrong!

% correct form: f = D Del2(x) + Grad(D) . Grad(x)
Del2x = Calc_Del2(x, Param.dx,Param.dy);
Gradx = Calc_Grad(x,Param.dx,Param.dy,1);
GradD = Calc_Grad(D,Param.dx,Param.dy,1);
GradD_dot_Gradx = Gradx(:,:,1) .* GradD(:,:,1) + Gradx(:,:,2) .* GradD(:,:,2);

f = D .* Del2x + GradD_dot_Gradx ;

