% OP_CURLU_CURLV_TP: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon curl u_j, curl v_i), exploiting the tensor product structure.
%
%   mat = op_curlu_curlv_tp (spu, spv, msh, [epsilon]);
%   [rows, cols, values] = op_curlu_curlv_tp (spu, spv, msh, [epsilon]);
%
% The same function works for 2d (scalar-valued curl) and 3d problems (vector-valued curl).
%
% INPUT:
%
%   spu:     object that defines the space of trial functions (see sp_vector)
%   spv:     object that defines the space of test functions (see sp_vector)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_cartesian)
%   epsilon: function handle to compute some physical coefficient (optional)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
%
% editor: Anna Ziegler

function varargout = op_dtC_dtA_curlu_curlv_tp(space1, space2, msh, msh_old, geo_map_jac_dCdt, t, coeff, matrix)

if isnan(geo_map_jac_dCdt)
    clear geo_map_jac_dCdt
    for idx = 1:msh.nel_dir(1)
        geo_map_jac_dCdt(:, :, :, :, idx) = NaN;
    end
end
for icomp = 1:space1.ncomp_param
    for idim = 1:msh.ndim
        size1 = size(space1.scalar_spaces{icomp}.sp_univ(idim).connectivity);
        size2 = size(space2.scalar_spaces{icomp}.sp_univ(idim).connectivity);
        if (size1(2) ~= size2(2) || size1(2) ~= msh.nel_dir(idim))
            error('One of the discrete spaces is not associated to the mesh')
        end
    end
end

A = spalloc(space2.ndof, space1.ndof, 3*space1.ndof);

for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col(msh, iel);
    msh_old_col = msh_evaluate_col(msh_old, iel);
    sp1_col = sp_evaluate_col(space1, msh_old_col, 'value', true, 'curl', true);
    sp2_col = sp_evaluate_col(space2, msh_old_col, 'value', true, 'curl', true);
    %     sp1_col = sp_evaluate_col_param (space1, msh_old_col, 'value', true, 'curl', true);
    %     sp2_col = sp_evaluate_col_param (space2, msh_old_col, 'value', true, 'curl', true);

    if (nargin >= 6)
        for idim = 1:msh.rdim
            x{idim} = reshape(msh_col.geo_map(idim, :, :), msh_col.nqn, msh_col.nel);
        end
        coeffs = coeff(x{:});
    else
        coeffs = ones(msh_col.nqn, msh_col.nel);
    end

    if (strcmp(matrix, 'stiffness'))
        A = A + op_dtC_curlu_curlv(sp1_col, sp2_col, msh_col, msh_old_col, geo_map_jac_dCdt(:, :, :, :, iel), t, coeffs);
    elseif (strcmp(matrix, 'mass'))
        A = A + op_dtA_u_v(sp1_col, sp2_col, msh_col, msh_old_col, geo_map_jac_dCdt(:, :, :, :, iel), t, coeffs);
    else
        error('Unknown matrix type, choose stiffness or mass');
    end
end

if (nargout == 1)
    varargout{1} = A;
elseif (nargout == 3)
    [rows, cols, vals] = find(A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
end

end
