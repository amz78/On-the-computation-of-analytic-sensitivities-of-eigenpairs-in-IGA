% OP_DTA_U_V_3D: assemble the matrix A = [a(i,j)], a(i,j) = (coeff curl u_j, curl v_i), with vector-valued curl.
%
%
% INPUT:
%
%   spu:    structure representing the space of trial functions (see sp_vector/sp_evaluate_col)
%   spv:    structure representing the space of test functions  (see sp_vector/sp_evaluate_col)
%   msh:    structure containing the domain partition and the quadrature
%               rule (see msh_cartesian/msh_evaluate_col) of the deformed mesh
%   msh_old:structure containing the domain partition and the quadrature
%               rule (see msh_cartesian/msh_evaluate_col) of the deformed mesh
%   jac_geo:either provided or if it is NaN, it will be computed
%   t:      deformation parameter in [0,1]
%   coeff:  physical parameter
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
%
% editor: Anna Ziegler

% this function is only implemented to work on the initial geometric
% domain! not on the parametric domain!

function varargout = op_dtA_u_v(spu, spv, msh, msh_old, jac_geo_iel, t, coeff)

global DERC
global DERA
global NI
global NMAX

rows = zeros(msh.nel*spu.nsh_max*spv.nsh_max, 1);
cols = zeros(msh.nel*spu.nsh_max*spv.nsh_max, 1);
values = zeros(msh.nel*spu.nsh_max*spv.nsh_max, 1);

quad_weights = msh_old.jacdet .* msh_old.quad_weights .* coeff;

ncounter = 0;
for iel = 1:msh.nel
    if (all(msh.jacdet(:, iel)))
        curlu_iel = reshape(spu.shape_functions(:, :, :, iel), 3, msh.nqn, 1, spu.nsh_max);
        curlv_iel = reshape(spv.shape_functions(:, :, :, iel), 3, msh.nqn, spv.nsh_max, 1);
        jac_geo_iel = zeros(size(msh.geo_map_jac(:, :, :, iel)-msh_old.geo_map_jac(:, :, :, iel))); %dV = Jacobian of old mesh - Jacobian of current mesh

        for iqn = 1:msh.nqn
            jac_geo_iel(:, :, iqn) = (msh.geo_map_jac(:, :, iqn, iel) - msh_old.geo_map_jac(:, :, iqn, iel)) * inv(msh_old.geo_map_jac(:, :, iqn, iel));
            dG_iel(:, :, iqn) = eye(size(jac_geo_iel(:, :, iqn))) + t * jac_geo_iel(:, :, iqn); %dG = I + t*dV
        end
        jacdet_iel = reshape(quad_weights(:, iel), [1, msh.nqn, 1, 1]);

        %memory pre-allocation
        trdVdG = zeros(msh.nqn, 1);
        detdG_dGT = zeros(size(dG_iel));
        dG_dV = zeros(size(jac_geo_iel));
        A = zeros(size(dG_iel));
        jacdet_curlu = bsxfun(@times, jacdet_iel, curlu_iel);

        for iqn = 1:msh.nqn
            detdG_dGT(:, :, iqn) = det(dG_iel(:, :, iqn)) * inv(dG_iel(:, :, iqn)); %det(dG) * dG^-1
            dG_dV(:, :, iqn) = (dG_iel(:, :, iqn)) \ jac_geo_iel(:, :, iqn); %(dG)^-1 * dV
            trdVdG(iqn) = trace(jac_geo_iel(:, :, iqn)/dG_iel(:, :, iqn)); %tr(dV dG^-1)

            A(:, :, iqn) = detdG_dGT(:, :, iqn) / dG_iel(:, :, iqn)'; %A = det(dG) * dG^-1 * dG^-T
            dG_dV_A(:, :, iqn) = dG_dV(:, :, iqn) * A(:, :, iqn); %dG^-1 * dV * A
            trdVdG_A(:, :, iqn) = trdVdG(iqn) * A(:, :, iqn); %tr(dV dG^-1) * A
        end

        tmp = zeros(msh.nqn, spu.nsh_max, spv.nsh_max);

        for iqn = 1:msh.nqn
            if isempty(DERA)
                dAdt(iqn, :, :) = bsxfun(@minus, trdVdG_A(:, :, iqn), bsxfun(@plus, dG_dV_A(:, :, iqn), dG_dV_A(:, :, iqn)')); %dAdt = -tr...
            else
                dnAdtn = DERA{NI};
                dViqn = jac_geo_iel(:, :, iqn);
                dAdt(iqn, :, :) = dnAdtn(dViqn(1, 1), dViqn(1, 2), dViqn(1, 3), dViqn(2, 1), dViqn(2, 2), dViqn(2, 3), dViqn(3, 1), dViqn(3, 2), dViqn(3, 3), t);
            end

            for nsh = 1:spu.nsh_max
                dAdT_curlu(:, iqn, :, nsh) = squeeze(dAdt(iqn, :, :)) * jacdet_curlu(:, iqn, :, nsh);
            end

            tmp(iqn, :, :) = sum(bsxfun(@times, dAdT_curlu(:, iqn, :, :), curlv_iel(:, iqn, :, :)), 1);
        end
        elementary_values = reshape(sum(tmp, 1), spv.nsh_max, spu.nsh_max);

        [rows_loc, cols_loc] = ndgrid(spv.connectivity(:, iel), spu.connectivity(:, iel));
        indices = rows_loc & cols_loc;
        rows(ncounter+(1:spu.nsh(iel) * spv.nsh(iel))) = rows_loc(indices);
        cols(ncounter+(1:spu.nsh(iel) * spv.nsh(iel))) = cols_loc(indices);
        values(ncounter+(1:spu.nsh(iel) * spv.nsh(iel))) = elementary_values(indices);
        ncounter = ncounter + spu.nsh(iel) * spv.nsh(iel);
    else
        warning('geopdes:jacdet_zero_at_quad_node', 'op_curlu_curlv_3d: singular map in element number %d', iel)
    end
end

if (nargout == 1)
    varargout{1} = sparse(rows(1:ncounter), cols(1:ncounter), ...
        values(1:ncounter), spv.ndof, spu.ndof);
elseif (nargout == 3 || nargout == 0)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
else
    error('op_curlu_curlv_3d: wrong number of output arguments')
end

end