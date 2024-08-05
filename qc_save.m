
function qc_save(parm, kind, J, lb, ub, zeta, u, u_analytic, switching, hamiltonian, varargin)

root = '../njp_aroqa/dat';  
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'path'
            root = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

data_filename = sprintf('%s/data_%s_norm%s_zeta%.0e.txt', root, parm, kind, zeta);
meta_filename = sprintf('%s/meta_%s_norm%s_zeta%.0e.txt', root, parm, kind, zeta);

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

