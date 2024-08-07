
x = 6; S = 208; N = 208; n = 100; zeta = 0.3;

errors = {};
for i = 1:4
    errors{i} = min(1, max(-1, randn(10, 1)));
end
        


data = [];
for i = S:N

    parm = sprintf('x%i_idx%.3i', x, i); u = {}; kinds = {'none', '2'};

    if exist(sprintf('dat/%s.csv', parm), 'file') == 2
        L = readmatrix(sprintf('dat/%s.csv', parm), 'Delimiter', ' ');
        [B, C] = problem_hamiltonian(L); F = B - C;
    else
        L = randn(x,x); [B, C] = problem_hamiltonian(L); F = B - C;
        writematrix(L, sprintf('dat/%s.csv', parm), 'Delimiter', ' ');
    end 


    
    [V,D] = eig(B);
    [~,I] = min(diag(D));
    x0    = V(:,I);
    u0    = 0.5*ones(n,1);


    for k = 1:size(kinds,2)

        ident = sprintf('%s_norm%s_zeta%.0e', parm, kinds{k}, zeta);
        if exist(sprintf('dat/data_%s.txt', ident), 'file') == 2
            [~,~,~,~,u{k}, ~,~,~,~] = qc_load(ident, 'path', 'dat');
        %else
            u0 = u{k};
            switch kinds{k}
                case {'none'}
                    [u{k}, ~, ~] = qc_optimize(C, F, u0, x0, kinds{k}, 'Zeta', zeta, ...
                        'Iterations', 300, 'Gamma', 0.2, 'Alpha', 0.8);
                    %[u{k}, ~, ~] = qc_optimize(C, F, u{k}, x0, kinds{k}, 'Zeta', zeta, ...
                    %    'Iterations', 700, 'Gamma', 0.4, 'Alpha', 0.3);
                case {'2', 'F', 'sqF', '2sh'}
                    [u{k}, ~, ~] = qc_optimize(C, F, u0, x0, kinds{k}, 'Zeta', zeta, ...
                        'Iterations', 300, 'Gamma', 0.2, 'Alpha', 0.8);
                otherwise
                    error('Invalid norm type supplied');
            end
            
        end

       
        [J, hamiltonian, u_analytic, switching, lb, ub] = qc_analyse(u{k}, C, F, x0, zeta, kinds{k});
        qc_save(parm, kinds{k}, J, lb, ub, zeta, u{k}, u_analytic, switching, hamiltonian, 'path', 'dat');
    end

    if exist(sprintf('dat/rob_%s.csv', parm), 'file') == 2
        II = readmatrix(sprintf('dat/rob_%s.csv', parm), 'Delimiter', ' ');
        Jm = II(:,2:3);
        Fm = II(:,3:5);
    else 
        [Bm,Jm,Fm] = eval_robustness(C, F, x0, u, kinds, errors);
        writematrix([Bm, Jm, Fm], sprintf('dat/rob_%s.csv', parm), 'Delimiter', ' ');
    end

    data = [data, Fm];

    hold on
    plot(1:size(Jm,1), Jm(:,1), 'r', 1:size(Jm,1), Jm(:,2), 'b');
    hold off

end

writematrix(data, sprintf('dat/x%i_result.csv', x), 'Delimiter', ' ');

function n_iter = select_iter(kind)
    switch kind
        case {'none'}
            n_iter = 1000;
        case {'2', 'F', 'sqF', '2sh'}
            n_iter = 300;
        otherwise
            error('Invalid norm type supplied');
    end
end

function gg = select_gamma(kind)
    switch kind
        case {'none'}
            gg = 0.8;
        case {'2', 'F', 'sqF', '2sh'}
            gg = 0.6;
        otherwise
            error('Invalid norm type supplied');
    end
end


function [Bm, Jm, Fm] = eval_robustness(C, F, x0, inputs, kinds, errors)

    Jm = [];
    Fm = [];
    Bm = [];

    for bb = 0.1:0.3:1.0
        er = {};
        for i = 1:size(errors,2)
            er{i} = bb*errors{i};
        end
        
        [Js, Fs, Bv] = qc_samples(C, F, x0, inputs, kinds, 'Error', er, 'MaxStep', 3e-2);
       
        Jmeans = [];
        Fmins  = [];
        for i = 1:size(kinds,2)
            Jmeans(i) = mean(Js{i});
            Fmins(i)  = min(Fs{i});
        end
    
        Jm = [Jm; Jmeans];
        Fm = [Fm; Fmins];
        Bm = [Bm; bb];
    end
   
end