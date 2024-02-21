function [sigma, subdif] = qc_norm(X, D, kind)
    switch kind 
        case {'none', 'qaoa'}
            sigma = 0;
            subdif = [0,0];
            
        case 'sqF'
            sigma = trace(X'*X);
            grad  = trace(D'*X + X'*D);
            subdif = [grad, grad];
            
        case 'F'
            sigma = sqrt(trace(X'*X));
            grad  = 0.5*trace(D'*X + X'*D)/sigma;
            subdif = [grad, grad];
       
            
        case '2shfd'
            epsi = 0.001;
            
            r1 = eig(X+epsi*D);
            o1 = -(max(r1)+min(r1))/2;
            r2 = eig(X-epsi*D);
            o2 = -(max(r2)+min(r2))/2;
 
            [~,S1,~] = svd(X+o1*eye(size(X,1))+epsi*D);
            [~,S2,~] = svd(X+o2*eye(size(X,1))-epsi*D);
            
            s1 = max(diag(S1));
            s2 = max(diag(S2));
            
            grad = (s1 - s2)/(2*epsi);
            sigma = 0.5* (s1 + s2);
            subdif = [grad, grad];

        case '2sh'
           
            [V,S] = eig(X);
            
            sigma = 0.5*( max(diag(S)) - min(diag(S)) );
            
            masks_1 = abs(diag(S)-max(diag(S))) < 1e-5;
            masks_2 = abs(diag(S)-min(diag(S))) < 1e-5;
            
            V_1 = V(:,masks_1);
            V_2 = V(:,masks_2);
            
            A_1 = 0.5*diag( V_1'*D*V_1 );
            A_2 = 0.5*diag( V_2'*D*V_2 );
            
            
            subdif = [min(real(A_1))-max(real(A_2)), max(real(A_1)) - min(real(A_2))];

            
        case '2'
            [U,S,V] = svd(X);
           
            sigma = max(diag(S));
            masks = abs(diag(S)-sigma) < 1e-5;
    
            U = U(:,masks);
            V = V(:,masks);
    
            A = diag( U'*D*V );

            subdif = [min(real(A)), max(real(A))];
            
        otherwise
            error('Invalid norm type supplied');
    end
end