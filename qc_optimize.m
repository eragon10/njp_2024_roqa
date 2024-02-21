    function [input, J, xe] = qc_optimize(C, F, u0, x0, kind, varargin) %(C, F, t0, x0, varargin) 

    m = size(C, 1);
    N = size(u0,1);
    %H = B - C;
    
    zeta  = 1e-2;
    niter = 100;
    gamma = 0.01;
    alpha = 0.0;
    T     = 1;
    
    
    
    %%
    spc_number = 1000;
    spc_values = zeros(spc_number+1,1);
    spc_sdiff  = zeros(spc_number+1,1);
    spc_index  = zeros(spc_number+1,1);
    norm_comp  = @(uin) nrm( uin );
            
    switch kind
        case {'qaoa'}
            [input, J, xe] = qc_qaoa(C, F, u0, x0, varargin{:});
            return;
        case {'2', '2sh', '2shfd'}
            for i = 1:(spc_number+1)
                [sig, dif] = qc_norm(C + ((i-1)/spc_number) * F, F, kind);
                spc_values(i) = sig;
                spc_sdiff(i)  = 0.5*sum(dif);
                spc_index(i)  = (i-1)/spc_number;
            end
            norm_comp = @( uin ) spc( uin );
        case {'none', 'F', 'sqF'}
        otherwise
            error('Invalid norm type supplied');
    end
    %%
    
    
    while ~isempty(varargin)
        switch lower(varargin{1})
            case 'time'
                T = varargin{2};
            case 'iterations'
                niter = varargin{2};
            case 'gamma'
                gamma = varargin{2};
            case 'alpha'
                alpha = varargin{2};
            case 'zeta'
                zeta = varargin{2};
            otherwise
                error(['Unexpected option: ' varargin{1}])
         end
         varargin(1:2) = [];
    end
    
    dt = T / N;
    
    
    input = u0;
    im = 0 .* u0;
    
    for k = 1:niter
        
        [~, dJdu, ~] = obj(input);
        
        if mod(k-1,10) == 0
           fprintf('Optimality condition: %0.5e\n', norm( (max(0, min(1, input - 0.01*dJdu) ) - input)/0.01 )  );
           plot(linspace(0,1,N), input);
           drawnow
        end
        
        
        im = alpha * im - gamma * real(dJdu);
        
        input = input + im;
        input = max(0, min(1, input) );
    end
    
    [J, ~, xe] = obj(input);
    
    function [J, dJdu, xe] = obj(u)
       dJdx = zeros(m,N);
       dJdu = zeros(N,1); 
       
       
       xe = x0;
       J  = 0;
       for j = 1:N
           Z = C + u(j) * F; 
           
           %E = expm(-1i*Z*dt);
           %D = -1i*F*E;
           E = eye(m)  - 1i * Z * dt - 1/2 * Z^2 * dt^2 + 1i/6 * Z^3 * dt^3;
           D = 1i/6*(F*Z^2 + Z*F*Z + Z^2*F)*dt^3 - 1/2*(F*Z + Z*F)*dt^2 - 1i * F * dt;
           

           %[V, dVdh] = fro(Z, H);
           %[V, dVdh] = spc( u(j) );
           [V, dVdh] = norm_comp( u(j) );
                 
           dJdx(:,1:(j-1)) = E * dJdx(:,1:(j-1));
           dJdx(:,j) = D * xe;
           
           dJdu(j)   = dVdh * dt * zeta;
           
           J  = J + V * dt * zeta;
           xe = E * xe;
       end
       
       dJdu = dJdu + (dJdx' * (C * xe)) + (dJdx.' * (C.' * conj(xe)));
       %dJdu = real( dJdu );
       
       
       J = J + real( xe' * C * xe );
    end

    function [V, dVdu] = spc( uin )
        V = interp1(spc_index, spc_values, uin ,'linear');
        dVdu = interp1(spc_index, spc_sdiff, uin ,'linear'); 
    end

    function [V, dVdu] = nrm( uin )
        [sigm, sdif] = qc_norm(C + uin * F, F, kind); 
        
        V = sigm;
        dVdu = 0.5*sum(sdif);
    end

end