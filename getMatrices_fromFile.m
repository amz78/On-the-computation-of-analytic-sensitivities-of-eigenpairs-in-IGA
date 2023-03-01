function [A, B, int_dofs, drchlt_dofs, geo_cavity, sp_cavity, msh_cavity] = getMatrices_fromFile(geo_cavity_file, deg, sub, dispDetailed)

[geo_cavity, bnd_cavity, intrfc_cavity, ~, boundary_interfaces] = mp_geo_load(geo_cavity_file);
npatch = numel(geo_cavity);

degree = [deg, deg, deg]; % Degree of the bsplines
regularity = [deg - 1, deg - 1, deg - 1]; % Regularity of the splines
nsub = [sub, sub, sub]; % Number of subdivisions
nquad = degree + 1; % Points for the Gaussian quadrature rule

% Physical parameters
eps0 = 8.85418781762e-12;
mu0 = 16 * atan(1) * 1.e-7;
% c     = 1 / sqrt (mu0 * eps0);
c_elec_perm = @(x, y, z) eps0 * ones(size(x));
c_magn_perm = @(x, y, z) mu0 * ones(size(x));

% -------------------------------------------------------------------------
% Construct the mesh and the space
% -------------------------------------------------------------------------

if dispDetailed
    fprintf('Building the mesh and the space:\n');
end

msh_cavity = cell(1, npatch);
sp = cell(1, npatch);
for iptc = 1:npatch
    [knots, zeta] = kntrefine(geo_cavity(iptc).nurbs.knots, nsub-1, degree, regularity);
    [knots_hcurl, degree_hcurl] = knt_derham(knots, degree, 'Hcurl');

    % Construct msh structure
    rule = msh_gauss_nodes(nquad);
    [qn, qw] = msh_set_quad_nodes(zeta, rule);
    msh_cavity{iptc} = msh_cartesian(zeta, qn, qw, geo_cavity(iptc));

    % Construct space structure
    scalar_spaces = cell(msh_cavity{iptc}.ndim, 1);
    for idim = 1:msh_cavity{iptc}.ndim
        scalar_spaces{idim} = sp_bspline(knots_hcurl{idim}, degree_hcurl{idim}, msh_cavity{iptc});
    end
    sp{iptc} = sp_vector(scalar_spaces, msh_cavity{iptc}, 'curl-preserving');
    clear scalar_spaces
end

msh_cavity = msh_multipatch(msh_cavity, bnd_cavity);
sp_cavity = sp_multipatch(sp, msh_cavity, intrfc_cavity, boundary_interfaces);
clear sp

% Compute total number of elements
nel_cavity = msh_cavity.nel;
ndof_cavity = sp_cavity.ndof;

if dispDetailed
    fprintf('\tnel : %i\n', nel_cavity);
    fprintf('\tndof: %i\n', ndof_cavity);
end

% -------------------------------------------------------------------------
% Assemble the matrices and apply boundary conditions
% -------------------------------------------------------------------------
if dispDetailed
    fprintf('Assemble the matrices and apply BC... ');
end

invmu = @(x, y, z) 1 ./ c_magn_perm(x, y, z);
A = op_curlu_curlv_mp(sp_cavity, sp_cavity, msh_cavity, invmu);
B = op_u_v_mp(sp_cavity, sp_cavity, msh_cavity, c_elec_perm);

% Apply homogeneous Dirichlet boundary conditions
Nbnd = cumsum([0, bnd_cavity.nsides]);
boundary_gnum = sp_cavity.boundary.gnum;
bnd_dofs = [];
for iref = 1:numel(bnd_cavity)
    iref_patch_list = Nbnd(iref) + 1:Nbnd(iref+1);
    bnd_dofs = union(bnd_dofs, [boundary_gnum{iref_patch_list}]);
end
drchlt_dofs = sp_cavity.boundary.dofs(bnd_dofs);
int_dofs = setdiff(1:ndof_cavity, drchlt_dofs);

if dispDetailed
    fprintf(' done\n');
end

%Reduce to DoFs
A = A(int_dofs, int_dofs);
B = B(int_dofs, int_dofs);

%Apply Scaling (Scaling factor s=1e10):
A = A ./ 1e10;
B = B * 1e10;
end
