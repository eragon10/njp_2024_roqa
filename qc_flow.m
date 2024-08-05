function [t, state, J] = qc_flow(C, F, x0, u, kind, varargin)
    
    T = 1;
    e = 0;
    maxstep = 1e-2;
    
    while ~isempty(varargin)
        switch lower(varargin{1})
            case 'time'
                T = varargin{2};
            case 'error'
                e = varargin{2};
            case 'maxstep'
                maxstep = varargin{2};
            otherwise
                error(['Unexpected option: ' varargin{1}])
        end
        varargin(1:2) = [];
    end
    
    
    cfg = odeset('RelTol',1e-6,'AbsTol',1e-10, 'MaxStep', maxstep);

    if size(e,1) > 1
        ej = @(t) interp1(linspace(0,T,size(e,1)), e, t, 'previous');
    else
        ej = @(t) e;
    end
     
    switch kind
        case {'2', '2sh', '2shfd', 'sqF', 'F', 'none'}
            uj = @(t) interp1(linspace(0,T,size(u,1)), u, t);
        case 'qaoa'
            uj = @(t) mod(sum(cumsum(u) < t), 2);
        otherwise
            error('Invalid kind');
    end
        
    switch kind
        case {'2sh', '2shfd'}
            phi = zeros(size(u,1),1);
            for qq = 1:size(phi,1)
                eqq = eig(C+F*u(qq));
                phi(qq) = -(max(eqq)+min(eqq))/2;
            end
                
            % + phij(t)*eye(size(C,1))
            phij = @(t) interp1(linspace(0,T,size(phi,1)), phi, t);
            fx = @(t,x) -1i * (1 + ej(t)) * (C + uj(t) * F + phij(t)*eye(size(C,1)) ) * x;
        otherwise
            fx = @(t,x) -1i * (1 + ej(t)) * (C + uj(t) * F) * x;
    end
        
    [t,y] = ode45(@(t,x) fx(t,x), [0,T], x0, cfg);
    
    J = (y') .* (C*y.');
    J = sum(real(J), 1).';
    state = y(end, :).';
    %J     = real(state'*C*state);

end