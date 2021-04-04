function Phi = ApplyBC_AC(Phi)

Phi(1,:) = 1;
Phi(:,1) = 1;
Phi(:,end) = 1;
Phi(end,:) = 1;

