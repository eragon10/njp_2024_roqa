function [Js, Fs, Bv] = qc_samples(C, F, x0, u0, opts, varargin)

T = 1;
E = { 0.1, -0.1 };

 maxstep = 1e-2;

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'time'
            T = varargin{2};
        case 'error'
            E = varargin{2};
        case 'maxstep'
            maxstep = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end


assert(size(u0, 2) == size(opts, 2), 'Should be the same size');

cfg = odeset('RelTol',1e-6,'AbsTol',1e-10, 'MaxStep', maxstep);
        
Js = cell(size(u0,2),1);
Fs = cell(size(u0,2),1);

Bs = cell(size(u0,2),1);
Bv = cell(size(u0,2),1);

for i = 1:size(u0, 2)
    [Bs{i}, Bv{i}] = rollout(u0{i}, 0, opts{i});
end



for k = 1:size(E,2)
    for i = 1:size(u0, 2)
        [xT, J] = rollout(u0{i}, E{k}, opts{i});
        fidelity = abs(xT'*Bs{i});
        Js{i} = [Js{i}, (J-Bv{i})];
        Fs{i} = [Fs{i}, fidelity];
    end
end


    function [state, J] = rollout(u,e, opt)        
        if size(e,1) > 1
            ej = @(t) interp1(linspace(0,T,size(e,1)), e, t, 'previous');
        else
            ej = @(t) e;
        end
        
        switch opt
            case {'2', '2sh', '2shfd', 'sqF', 'F', 'none'}
                uj = @(t) interp1(linspace(0,T,size(u,1)), u, t);
            case 'qaoa'
                uj = @(t) mod(sum(cumsum(u) < t), 2);
            otherwise
                error('Invalid kind');
        end
        
        switch opt
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
        
        [~,y] = ode45(@(t,x) fx(t,x), [0,T], x0, cfg);
        state = y(end, :).';
        J     = real(state'*C*state);
    end

end