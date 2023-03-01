function [dA, dB] = getDerivativeMatrices(space_pb, mesh_pb, mesh_tesla, int_dofs_pb, t)
%computes derivatives of system matrices up to NMAXth derivative

% INPUT:
%
% space_pb: space of the undeformed geometry
% mesh_pb:  mesh of the undeformed geometry
% mesh_tesla:  mesh of the deformed geometry
% int_dofs_pb: internal DoFs
% t:        parameter at which to compute the derivatives
%
% author: Anna Ziegler
% E-Mail: anna.ziegler@tu-darmstadt.de

global NI
global NMAX

eps0 = 8.85418781762e-12;
mu0 = 16 * atan(1) * 1.e-7;
c_elec_perm = @(x, y, z) eps0 * ones(size(x));
c_magn_perm = @(x, y, z) mu0 * ones(size(x));
invmu = @(x, y, z) 1 ./ c_magn_perm(x, y, z);

dA = cell(NMAX, 1);
dB = cell(NMAX, 1);

%%
jac_V = NaN;

% compute derivatives of system matrices
for NI = 1:NMAX
    fprintf('Compute Derivative %i of Stiffness Matrix \n', NI);
    dK_tesla = op_dtC_dtA_curlu_curlv_mp(space_pb, space_pb, mesh_tesla, mesh_pb, jac_V, t, invmu, 'stiffness');
    dK_tesla = dK_tesla(int_dofs_pb, int_dofs_pb);
    dA{NI} = dK_tesla ./ 1e10;
    clear dK_tesla

    fprintf('Compute Derivative %i of Mass Matrix \n', NI);
    dM_tesla = op_dtC_dtA_curlu_curlv_mp(space_pb, space_pb, mesh_tesla, mesh_pb, jac_V, t, c_elec_perm, 'mass');
    dM_tesla = dM_tesla(int_dofs_pb, int_dofs_pb);
    dB{NI} = dM_tesla .* 1e10;
    clear dM_tesla
end

NI = 1;
end
