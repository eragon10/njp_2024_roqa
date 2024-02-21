function [input,J,xe] = qc_qaoa(C, F, t0, x0, varargin) 

    m = size(C, 1);
    N = size(t0,1);
    T = sum(t0);
    
    niter = 100;
    gamma = 0.005;
    alpha = 0.0;
    
    while ~isempty(varargin)
        switch lower(varargin{1})
            case 'zeta'
            case 'time'
                T = varargin{2};
            case 'iterations'
                niter = varargin{2};
            case 'gamma'
                gamma = varargin{2};
            case 'alpha'
                alpha = varargin{2};
            otherwise
                error(['Unexpected option: ' varargin{1}])
         end
         varargin(1:2) = [];
    end
    
    [Vc, Dc] = eig(C, 'vector');
    [Vb, Db] = eig(C + F, 'vector');
    
    input = t0;
    veloc = 0 .* t0;
    
    for k = 1:niter
   
        [~, dJdu, ~] = obj(input);
        if mod(k+9, 10) == 0
            fprintf('Optimality condition: %0.5e\n', ...
                norm(simplex_projection(input - 0.01*dJdu, 'Iterations', 500, 'Scale', T) - input)/0.01);
           
        
            xx = []; % zeros(2*N,1);
            yy = []; %zeros(2*N,1); 
            ts = 0;
            for i = 1:N
                %xx(2*i-1) = ts;
                %yy(2*i-1) = mod(i+1,2);
                if abs(input(i)) < 1e-6
                    continue;
                end
                
                xx = [xx, ts, ts + input(i)];
                yy = [yy, mod(i+1,2), mod(i+1,2)];
                ts = ts + input(i);
               
                %xx(2*i)   = ts;
                %yy(2*i)   = mod(i+1,2);
            end
            %yy = mod( sum(tril(ones(N,N)) * input >= xx, 1), 2 );
            plot(xx, yy);
            xlim([-0.1,1.1]);
            ylim([-0.1,1.1]);
            drawnow;
        end
       
        %dJdu
        
        veloc = alpha * veloc - gamma * real(dJdu);
        
        input = input + veloc;
        
        input = simplex_projection(input, 'Iterations', 500, 'Scale', T);
        
    end
   
    [J, ~, xe] = obj(input);
    
    function [J, dJdt, xe] = obj(t)
       dPdt = zeros(m,N);
       
       xe = x0;
       for j = 1:N
           switch mod(j,2)  
               case 0 
                E = Vb * diag(exp(-1i*Db .* t(j))) * Vb';
                D = -1i * Vb * diag(Db .* exp(-1i*Db .* t(j))) * Vb';
               case 1
                E = Vc * diag(exp(-1i*Dc .* t(j))) * Vc';
                D = -1i * Vc * diag(Dc .* exp(-1i*Dc .* t(j))) * Vc';     
           end
           
                
           dPdt(:,1:(j-1)) = E * dPdt(:,1:(j-1));       
           dPdt(:,j) = D * xe;
           
           xe = E * xe;
           xe = xe / sqrt(xe'*xe);  
       end
       
       %dPdt
       J = real( xe' * C * xe );
       
       dJdt = (dPdt' * (C * xe)) + (dPdt.' * (C.' * conj(xe)));
       dJdt = real( dJdt );
    end

end