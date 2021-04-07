function A = ApplyBC(A)

w = 10;
T0 = 273+2;
A(1:w,:) = T0;
A(end-w:end,:) = T0;

A(:,1:w) = T0;
A(:,end-w:end) = T0;


end