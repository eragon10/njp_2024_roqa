function [B,C,F,x0] = x4params()

    L = [    0.8147    0.6324    0.9575    0.9572;
             0.9058    0.0975    0.9649    0.4854;
             0.1270    0.2785    0.1576    0.8003;
             0.9134    0.5469    0.9706    0.1419 ];
    
    [B, C] = problem_hamiltonian(L);
    
    F = B - C;
    
    
    [V,D] = eig(B);
    [~,I] = min(diag(D));
    x0    = V(:,I);


end