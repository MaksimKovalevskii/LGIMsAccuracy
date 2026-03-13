% SEEMS THAT WRONG FUNCTION!!!+Bs for global system, -Bs for local

function Ghat= GmatrixEP(x, y, z, b0, b1, b2, b3)
%syms x y z b0 b1 b2 b3

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
b = [b1 b2 b3];
Bs = skew(b);

Ghat = 2*[-b.' (b0*eye(3)-Bs)];
end