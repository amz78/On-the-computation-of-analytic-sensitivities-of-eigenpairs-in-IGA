function [eigf_res, eigv_res] = solveMyEVP(sys, neig, sigma)

%% function to solve generalized eigenvalueproblems for sys
%
% INPUT:
% sys:  structure 'sys': stiffness matrix sys.A; massmatrix sys.B
% and derivatives sys.dA and sys.dB
% neig: number of computed eigenvalues
% sigma: scalar value used as specified in eigs -> search eigenvalues
%        closest to sigma

%% OUTPUT
% eigv_res: vector of filtered and sorted eigenvalues
% eigf_res: matrix of filtered (same as eigv_res) and sortedt eigenvalues

%%


sys.A = 0.5 * (sys.A + sys.A');
sys.B = 0.5 * (sys.B + sys.B');

opts.maxit = 500;
opts.tol = 1e-8;

opts.spdB = true;
[eigv, eigf, flag] = eigs(sys.A, sys.B, neig, sigma/1e0, opts);


if flag
    warning('eigs did not converge. Please check results.');
end
[eigf, perm] = sort(diag(eigf));
eigv = eigv(:, perm);
n = numel(find(eigf < 1e15/1e20));
eigf_res = eigf(n+1:end);
eigv_res = eigv(:, n+1:end);
end
