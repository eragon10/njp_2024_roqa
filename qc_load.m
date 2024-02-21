

function [J, lb, ub, zeta, u, u_analytic, switching, hamiltonian, kind] = qc_load(filename)

kind = extractBetween(filename,'norm', '_');
kind = kind{:};

data = num2cell(readmatrix(sprintf('../dat/data_%s.txt', filename), 'Delimiter', ' '), 1);
meta = num2cell(readmatrix(sprintf('../dat/meta_%s.txt', filename), 'Delimiter', ' '), 1);

[J, lb, ub, zeta] = deal(meta{:});

switch kind
    case {'2', '2sh', 'F', 'none'}
        [u, u_analytic, switching, hamiltonian] = deal(data{:});
    case {'qaoa'}
        [u, ~] = deal(data{:});
        u = u(2:end);
        u_analytic  = 0.*u;
        switching   = 0.*u;
        hamiltonian = 0.*u;
    otherwise
        error('Invalid kind supplied');
end

end

