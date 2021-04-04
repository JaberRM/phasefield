function E = Calc_Strain(U,Param)

E = zeros(size(U,1) , size(U,2) ,3);

GradU1 = Calc_Grad(U(:,:,1),Param.dx,Param.dy,1);
GradU2 = Calc_Grad(U(:,:,2),Param.dx,Param.dy,1);


E(:,:,1) = GradU1(:,:,1);
E(:,:,2) = GradU2(:,:,2);
E(:,:,3) = 0.5 * ( GradU1(:,:,2) + GradU2(:,:,1) );




