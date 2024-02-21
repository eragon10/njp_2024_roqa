function [B,C,F,x0] = x3params()

    L = [      0.9298    0.9143   -0.7162;
              -0.6848   -0.0292   -0.1565;
               0.9412    0.6006    0.8315 ];
    
    [B, C] = problem_hamiltonian(L);
    
    F = B - C;
    
    
    [V,D] = eig(B);
    [~,I] = min(diag(D));
    x0    = V(:,I);


end