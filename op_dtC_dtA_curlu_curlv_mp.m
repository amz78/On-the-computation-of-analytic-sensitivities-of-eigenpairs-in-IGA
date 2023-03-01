% OP_DTC_CURLU_CURLV_MP: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon curl u_j, curl v_i), in a multipatch domain.
%
%   mat = op_dtC_curlu_curlv_mp (spu, spv, msh, [epsilon], [patches])
%
% The same function works for 2d (scalar-valued curl) and 3d problems (vector-valued curl).
%
% INPUT:
%
%   spu:     object that defines the space of trial functions (see sp_multipatch)
%   spv:     object that defines the space of test functions (see sp_multipatch)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_multipatch)
%   epsilon: function handle to compute some physical coefficient. Equal to one if left empty.
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%
% editor: Anna Ziegler

function A = op_dtC_dtA_curlu_curlv_mp(spu, spv, msh, msh_old, geo_map_jac_dCdt, t, coeff, matrix, patch_list)

if (nargin < 9)
    patch_list = 1:msh.npatch;
end

if (msh.npatch ~= msh_old.npatch)
    error('msh and msh_old do not have the same number of patches')
end

if ((spu.npatch ~= spv.npatch) || (spu.npatch ~= msh.npatch))
    error('op_curlu_curlv_mp: the number of patches does not coincide')
end

if ~iscell(geo_map_jac_dCdt)
    if isnan(geo_map_jac_dCdt)
        clear geo_map_jac_dCdt
        for idx = patch_list
            geo_map_jac_dCdt{idx} = NaN;
        end
    end
end

ncounter = 0;
for iptc = patch_list
    if (nargin < 7 || isempty(coeff))
        [rs, cs, vs] = op_dtC_dtA_curlu_curlv_tp(spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, msh_old.msh_patch{iptc}, geo_map_jac_dCdt{iptc}, t, matrix);
    else
        [rs, cs, vs] = op_dtC_dtA_curlu_curlv_tp(spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, msh_old.msh_patch{iptc}, geo_map_jac_dCdt{iptc}, t, coeff, matrix);
    end
    rows(ncounter+(1:numel(rs))) = spv.gnum{iptc}(rs);
    cols(ncounter+(1:numel(rs))) = spu.gnum{iptc}(cs);

    if (~isempty(spv.dofs_ornt))
        vs = spv.dofs_ornt{iptc}(rs)' .* vs;
    end
    if (~isempty(spu.dofs_ornt))
        vs = vs .* spu.dofs_ornt{iptc}(cs)';
    end

    vals(ncounter+(1:numel(rs))) = vs;
    ncounter = ncounter + numel(rs);
end

A = sparse(rows, cols, vals, spv.ndof, spu.ndof);
clear rows cols vals rs cs vs

end