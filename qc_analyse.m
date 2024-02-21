function [J, hamiltonian, u_analytic, switching, lb, ub] = qc_analyse(u, C, F, x0, zeta, kind)

N = size(u,1);
m = size(F,1);

dt = 1 / N;

cCF   = C*F - F*C;    
cCcCF = C * cCF - cCF * C;
cFcCF = F * cCF - cCF * F;


switch kind
    case {'qaoa'}
        lb = NaN;
        ub = NaN;
    case {'none'}
        lb = 0;
        ub = 0;
    case {'2', 'F', 'sqF', '2sh'}
        [~,dl] = qc_norm(C+0.01*F,F,kind);
        [~,du] = qc_norm(C+0.99*F,F,kind);
        lb = zeta*dl(2);
        ub = zeta*du(1);
    otherwise
        error('Invalid norm type supplied');
end

%Z     = C + u(j) * F;
%trans = eye(m)  - 1i * Z * dt - 1/2 * Z^2 * dt^2 + 1i/6 * Z^3 * dt^3;

hamiltonian = zeros(1,N);
u_analytic  = zeros(1,N);
switching   = zeros(1,N);

x = cell(N,1);
l = cell(N,1);
t = cell(N,1);

state = x0;

for i = 1:N
    switch kind
        case {'2', 'F', 'none'}
            t{i} = expm( -1i*(C+u(i)*F)*dt  );
        case {'2sh', '2shfd'}      
            eigen = eig(C+u(i)*F);
            off = -(max(eigen)+min(eigen))/2;
            %off = 0;
            t{i} = expm( -1i*(C+u(i)*F+off*eye(size(C,1)))*dt );
        case {'qaoa'}
            return;
        otherwise
            error('Invalid error supplied');
    end

    x{i} = t{i} * state;
    state = x{i} / sqrt(x{i}'*x{i});
end
J = real( state' * C * state );

% endpoint condition
state = - C * state;
scale = sqrt(state'*state);


switch kind
    case {'2', '2sh'} 
        nl = 400;
        xl = linspace(0,1,nl);
        yl = zeros(1,nl);

        for k = 1:nl
            [~, sub] = qc_norm(C+xl(k)*F,F, kind);
            yl(k) = 0.5*sum(sub);
        end
        spcnorm_subdif_inverse = @(x) interp1(yl,xl,x,'linear'); 
end

for i = N:-1:1
    
    hamiltonian(i) = switching_function(C+u(i)*F, x{i}, state) - zeta * qc_norm(C+u(i)*F,F,kind);
    switching(i)   = switching_function(F,x{i},state);
    
    switch kind
        case 'qaoa'
            u_analytic(i) = 0;
            hamiltonian(i) = 0;
            switching(i) = 0;
        case 'none'
            epsilon = 0.01;
            s = switching_function(F,x{i},state);
            if s > epsilon
                u_analytic(i) = 1;
            elseif s < -epsilon
                u_analytic(i) = 0;
            else
                u_analytic(i)  = min(1, max(0, real( ...
                    - switching_function(cCcCF, x{i}, state)/switching_function(cFcCF, x{i}, state) ...
                    )));
            end
        case 'sqF'
            u_analytic(i) = ( switching_function(F, x{i}, state)  ...
                - zeta * trace(C'*F + F'*C) ) / ( 2 * zeta * trace(F'*F) );
        case 'F'
            s = switching_function(F, x{i}, state) / zeta;
            if s >= ub/zeta
                u_analytic(i) = 1;
            elseif s <= lb/zeta
                u_analytic(i) = 0;
            else
                u_analytic(i) = s / sqrt(trace(F*F)) * sqrt( (trace(C*C)-trace(C*F)^2/trace(F*F)) / (trace(F*F) - s^2) ) ...
                    - trace(C*F) / trace(F*F);
            end
        case {'2', '2sh'}
            
            s = switching_function(F,x{i},state) / zeta;
            
            if s >= ub/zeta
                u_analytic(i) = 1;
            elseif s <= lb/zeta
                u_analytic(i) = 0;
            else
                u_analytic(i) = spcnorm_subdif_inverse( s );
            end
            %grad = switching_function(F,x{i},state);
            %u_analytic(i) = or(and(grad >= zeta * subdif(1), grad <= zeta*subdif(2)), ...
            %     or( abs( grad - zeta * subdif(1) ) < 1e-1, abs( grad - zeta * subdif(2) ) < 1e-1 ) ...
            %    ); 
         
    end
    
    
    l{i}  = t{i}'*state;
    state = scale * l{i} / sqrt(l{i}'*l{i});
    
end

    function v = switching_function(X, s, l)
        v = real( -1i*l'*X*s + 1i*s'*X*l );
    end

end