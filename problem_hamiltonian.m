function [B, C] = problem_hamiltonian(L)
    n = size(L,1);

    sigma_x = [0,1;1,0];
    %sigma_z = [0,-1;1,0];
    sigma_z = [1,0;0,-1];
    
    B = zeros(2^n, 2^n);
    for i = 1:n
        sig_i = [kron(ones(1,i-1),eye(2)), sigma_x, kron(ones(1,n-i), eye(2))];
        B = B - generate( sig_i );
    end
    
    C = zeros(2^n, 2^n);
    for i = 1:n
       for k = 1:n
           
           sig_i = [kron(ones(1,i-1),eye(2)), sigma_z, kron(ones(1,n-i), eye(2))];
           sig_k = [kron(ones(1,k-1),eye(2)), sigma_z, kron(ones(1,n-k), eye(2))];
           
           C = C + L(i,k) * generate(sig_i) * generate(sig_k);
       end
    end
    
    
    
    function [M] = generate(Z)
        switch n 
            case 1
                M = Z;
            otherwise
                M = kron(Z(1:2, 1:2), Z(1:2,3:4));
                for j = 3:n
                    M = kron(M, Z(1:2, (2*j-1):(2*j)));
                end
            
        end
    end
end