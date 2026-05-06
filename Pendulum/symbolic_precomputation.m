clear all;
clc;

syms x y z psi1 psi2 psi3 real
syms xd yd zd psi1d psi2d psi3d real

% Define generalised vector q, derivatives qd, rotational part p
q=[x y z psi1 psi2 psi3];
qd=[xd yd zd psi1d psi2d psi3d];
p = [psi1, psi2, psi3];

% Skew symmetric matrix
   function S = skew(v)
        S = [0 -v(3) v(2); 
             v(3) 0 -v(1); 
            -v(2) v(1) 0];
   end
Ps = skew (p);

%  Determining angle phi (theta in report) from psi
Phi = sqrt(psi1^2 + psi2^2 + psi3^2);
        
%  Transform matrix A and A2 (for small phi) from Rodriguez formula
A = eye(3) + Ps*(sin(Phi)/Phi) + 2*Ps*Ps*sin(Phi/2)*sin(Phi/2)/(Phi*Phi);
A2 = eye(3) + Ps*(1 - Phi^2/6 + Phi^4/120) + 2*Ps*Ps*(1/2 - Phi^2/48 + Phi^4/3840);

%  Constraints for spherical joint
rp = [0; 0; -3];
R = [x; y; z];
C = (R + A*rp).';
C2 = (R + A2*rp).';

% Symbolic Jacobian (Cq)
Cq = jacobian(C, q);
Cq2 = jacobian(C2, q);

% Symbolic Qc (forces) 
Cq_qd = Cq*qd';
Cq2_qd = Cq2*qd';

dCq_qd_dq = jacobian(Cq_qd, q);
dCq_qd_dq2 = jacobian(Cq2_qd, q);

Qc_sym = -dCq_qd_dq*qd';
Qc_sym2 = -dCq_qd_dq2*qd';

% simplifiying expressions due to lenghtiness 
Cq = simplify(Cq);
Cq2 = simplify(Cq2);
Qc_sym = simplify(Qc_sym);
Qc_sym2 = simplify(Qc_sym2);

% Generate MATLAB Functions - Jacobians, constraints C, dC (Cq_qd) for
% small phi and other phi
matlabFunction(Cq, 'File', 'Cq_func', 'Vars', {q}, 'Optimize', true);
matlabFunction(Cq2, 'File', 'Cq_zero', 'Vars', {q}, 'Optimize', true);
matlabFunction(Qc_sym, 'File', 'Qc_func', 'Vars', {q, qd}, 'Optimize', true);
matlabFunction(Qc_sym2, 'File', 'Qc_zero', 'Vars', {q, qd}, 'Optimize', true);

matlabFunction(C, 'File', 'C_func', 'Vars', {q}, 'Optimize', true);
matlabFunction(C2, 'File', 'C_zero', 'Vars', {q}, 'Optimize', true);
matlabFunction(Cq_qd, 'File', 'Cq_qd_func', 'Vars', {q, qd}, 'Optimize', true);
matlabFunction(Cq2_qd, 'File', 'Cq_qd_zero', 'Vars', {q, qd}, 'Optimize', true);

disp('Symbolic precomputation completed. Functions Cq_func.m,Cq_zero.m  and Qc_func.m, Qc_zero.m have been generated.');