function T = Calc_Stress(E,Param)

% T(:,:,1) = Param.C(:,:,1,1) * E(:,:,1) + Param.C(:,:,1,2) * E(:,:,2) + Param.C(:,:,1,3) * E(:,:,3);
% T(:,:,2) = Param.C(:,:,2,1) * E(:,:,1) + Param.C(:,:,2,2) * E(:,:,2) + Param.C(:,:,2,3) * E(:,:,3);
% T(:,:,3) = Param.C(:,:,3,1) * E(:,:,1) + Param.C(:,:,3,2) * E(:,:,2) + Param.C(:,:,3,3) * E(:,:,3);


T(:,:,1) = Param.C(1,1) * E(:,:,1) + Param.C(1,2) * E(:,:,2) + Param.C(1,3) * E(:,:,3);
T(:,:,2) = Param.C(2,1) * E(:,:,1) + Param.C(2,2) * E(:,:,2) + Param.C(2,3) * E(:,:,3);
T(:,:,3) = Param.C(3,1) * E(:,:,1) + Param.C(3,2) * E(:,:,2) + Param.C(3,3) * E(:,:,3);