
function qc_save(parm, kind, J, lb, ub, zeta, u, u_analytic, switching, hamiltonian)

data_filename = sprintf('../dat/data_%s_norm%s_zeta%.0e.txt', parm, kind, zeta);
meta_filename = sprintf('../dat/meta_%s_norm%s_zeta%.0e.txt', parm, kind, zeta);

writematrix([J, lb, ub, zeta], meta_filename, 'Delimiter', ' ');
        
switch kind
    case {'2', '2sh', 'F', 'none'}
        writematrix([u, u_analytic.', switching.', hamiltonian.'], data_filename, 'Delimiter', ' ');
    case {'qaoa'}
        u = [0; u];
        writematrix([u, tril(ones(size(u,1),size(u,1)))*u], data_filename, 'Delimiter', ' ');
    otherwise
        error('Invalid kind supplied');
end

end

