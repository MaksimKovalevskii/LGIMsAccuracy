function Ghat2= GmatrixCart(x, y, z, psi1, psi2, psi3)

% Skew symmetric matrix
function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

% rotational part of vector q ( p )
psi = [psi1 psi2 psi3];
PsiHat = skew(psi);

%  Determining angle phi (theta in report) from psi
phi=sqrt(psi1^2+psi2^2+psi3^2);

%  Matrix G from Rodriguez formula (with Taylor series for small phi)
if abs(phi) < 10^-3 
Ghat2 = eye(3) + PsiHat*(-0.5+phi^2/24-phi^4/720+phi^6/40320) + PsiHat*PsiHat*(1/6-phi^2/120+phi^4/5040-phi^6/362880);
else
Ghat2 = eye(3) - PsiHat*(1-cos(phi))/(phi^2) + PsiHat*PsiHat*(phi-sin(phi))/(phi^3);
end
end