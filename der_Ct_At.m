function [der] = der_Ct_At(derorder, fun)
%   Computes the derivatives of C(t) and A(t) symbolically
%   and returns the symbolic expressions of the derivatives up to
%   derorder-th order as matlab function handle

%   INPUT:
%   derorder :  order of highest derivative to compute
%   fun:        set to 'C' or 'A' to choose which derivatives to compute (C
%               for stiffness, A for mass matrix)
%
% author: Anna Ziegler
% E-Mail: anna.ziegler@tu-darmstadt.de

syms t dV11 dV12 dV13 dV21 dV22 dV23 dV31 dV32 dV33
dV = [dV11, dV12, dV13; dV21, dV22, dV23; dV31, dV32, dV33];
dG = eye(3) + t * dV;

if fun == 'C'
    C(t) = 1 / det(dG) * dG' * dG;
    fder(dV11, dV12, dV13, dV21, dV22, dV23, dV31, dV32, dV33, t) = C(t);

elseif fun == 'A'
    A(t) = det(dG) * inv(dG) * inv(dG');
    fder(dV11, dV12, dV13, dV21, dV22, dV23, dV31, dV32, dV33, t) = A(t);
else
    error('Invalid function to derive, must be ''C'' or ''A''. ')
end

der = cell(derorder, 1);

for i = 1:derorder
    fder(dV11, dV12, dV13, dV21, dV22, dV23, dV31, dV32, dV33, t) = diff(fder, t);
    der{i} = matlabFunction(fder);
end

end