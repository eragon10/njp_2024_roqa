function [B,C,F,x0] = x2params()

    L = [0.2 1;
         0.5 3.0];
    
    [B, C] = problem_hamiltonian(L);
    
   
    F = B - C;  
    F = 0.7*F;
    B = F + C;
    
    B = 0.5 * B;
    C = 0.5 * C;
    F = B - C;
    
    [V,D] = eig(B);
    [~,I] = min(diag(D));
    x0    = V(:,I);


end