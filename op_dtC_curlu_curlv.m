% OP_DTC_CURLU_CURLV_3D: assemble the matrix A = [a(i,j)], a(i,j) = (coeff curl u_j, curl v_i), with vector-valued curl.
%
%   mat = op_dtC_curlu_curlv_3d (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_curlu_curlv_3d (spu, spv, msh, epsilon);
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


function varargout = op_dtC_curlu_curlv(spu, spv, msh, msh_old, jac_geo, t, coeff)

global DERC
global DERA
global NI
global NMAX

%IF CHANGING THIS SETTING, ALSO ADAPT IN OP_DTC_DTA_CURLU_CURLV_TP.m !
%THERE: DON'T EVALUATE SPACE BY SP_EVALUATE_COL, BUT WITH
%SP_EVALUATE_COL_PARAM !
paramDomain = false; %choose on which domain to compute: either
%parametric or physical domain (msh old). Both work,
%but require different mappings

%check if a jac_geo to use exists or need to be computed
if isnan(jac_geo)
    compute_jac_geo = true;
    clear jac_geo;
else
    compute_jac_geo = false;
end
rows = zeros(msh.nel*spu.nsh_max*spv.nsh_max, 1);
cols = zeros(msh.nel*spu.nsh_max*spv.nsh_max, 1);
values = zeros(msh.nel*spu.nsh_max*spv.nsh_max, 1);

%set weights depending on domain
if paramDomain
    quad_weights = msh_old.quad_weights .* coeff;
else
    quad_weights = msh_old.jacdet .* msh_old.quad_weights .* coeff;
end

ncounter = 0;
for iel = 1:msh.nel
    if exist('jac_geo')
        jac_geo_iel = jac_geo(:, :, :, iel);
    end
    if (all(msh.jacdet(:, iel)))
        curlu_iel = reshape(spu.shape_function_curls(:, :, :, iel), 3, msh.nqn, 1, spu.nsh_max);
        curlv_iel = reshape(spv.shape_function_curls(:, :, :, iel), 3, msh.nqn, spv.nsh_max, 1);

        %if not using an existing jac_geo, it is computed here depending on the chosen domain
        if compute_jac_geo %corresponds to displacement vector field V
            for iqn = 1:msh.nqn
                if paramDomain %^x \in ^\Omega, then dV(^x) = dG2(^x) - dG1(^x)
                    jac_geo_iel(:, :, iqn) = (msh.geo_map_jac(:, :, iqn, iel) - msh_old.geo_map_jac(:, :, iqn, iel));
                else %apply inverse Jacobian for computation on msh_old, corresponds to: x \in \Omega_1, then dV(x) = (dG2(^x) - dG1(^x))*1/(dG1(^x))
                    jac_geo_iel(:, :, iqn) = (msh.geo_map_jac(:, :, iqn, iel) - msh_old.geo_map_jac(:, :, iqn, iel)) * inv(msh_old.geo_map_jac(:, :, iqn, iel));
                end
            end
        end
        for iqn = 1:msh.nqn
            %compute jacobian of mapping G (dG) on chosen domains
            if paramDomain
                dG_iel(:, :, iqn) = msh_old.geo_map_jac(:, :, iqn, iel) + t * jac_geo_iel(:, :, iqn); %dG(^x) = dG1(^x) + t*dV(^x)
            else
                dG_iel(:, :, iqn) = eye(size(jac_geo_iel(:, :, iqn))) + t * jac_geo_iel(:, :, iqn); %dG(x) = I + t*dV(x)
            end
        end
        jacdet_iel = reshape(quad_weights(:, iel), [1, msh.nqn, 1, 1]);

        %memory pre-allocation
        trdVdG = zeros(msh.nqn, 1);
        detdG_dGT = zeros(size(dG_iel));
        detdG_dVT = zeros(size(jac_geo_iel));
        C = zeros(size(dG_iel));
        jacdet_curlu = bsxfun(@times, jacdet_iel, curlu_iel);

        for iqn = 1:msh.nqn
            detdG_dGT(:, :, iqn) = 1 / det(dG_iel(:, :, iqn)) * dG_iel(:, :, iqn)'; %1/det(dG) * dG'
            detdG_dVT(:, :, iqn) = 1 / det(dG_iel(:, :, iqn)) * (jac_geo_iel(:, :, iqn))'; %1/det(dG) * dV'
            trdVdG(iqn) = trace((jac_geo_iel(:, :, iqn))/dG_iel(:, :, iqn)); %tr(dV dG^-1)

            C(:, :, iqn) = detdG_dGT(:, :, iqn) * dG_iel(:, :, iqn); %C = 1/det(dG) * dG' * dG
            invDet_dVT_dG(:, :, iqn) = detdG_dVT(:, :, iqn) * dG_iel(:, :, iqn); %1/det(dG) * (dV' * dG)
            trdVdG_C(:, :, iqn) = -trdVdG(iqn) * C(:, :, iqn); %-tr(dV dG^-1) * C
        end

        tmp = zeros(msh.nqn, spu.nsh_max, spv.nsh_max);

        for iqn = 1:msh.nqn
            if isempty(DERC)
                dCdt(iqn, :, :) = bsxfun(@plus, trdVdG_C(:, :, iqn), bsxfun(@plus, invDet_dVT_dG(:, :, iqn), invDet_dVT_dG(:, :, iqn)'));
            else

                dnCdtn = DERC{NI};
                dViqn = jac_geo_iel(:, :, iqn);
                newdCdt(iqn, :, :) = dnCdtn(dViqn(1, 1), dViqn(1, 2), dViqn(1, 3), dViqn(2, 1), dViqn(2, 2), dViqn(2, 3), dViqn(3, 1), dViqn(3, 2), dViqn(3, 3), t);

                dCdt(iqn, :, :) = newdCdt(iqn, :, :);
            end
            for nsh = 1:spu.nsh_max
                dCdT_curlu(:, iqn, :, nsh) = squeeze(dCdt(iqn, :, :)) * jacdet_curlu(:, iqn, :, nsh);
            end

            tmp(iqn, :, :) = sum(bsxfun(@times, dCdT_curlu(:, iqn, :, :), curlv_iel(:, iqn, :, :)), 1); %tmp = dCdt_curlNj_curlNi;
        end
        elementary_values = reshape(sum(tmp, 1), spv.nsh_max, spu.nsh_max); %"integrieren" über iqn

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