

n = 200;
m = 20;

zeta = 0.3;
kind = 'none';
parm = 'x6params_1';

%[B,C,F,x0] = x3params_extreme();
[B,C,F,x0] = eval(sprintf('%s()', parm));


[~,s1] = qc_norm(C+0.01*F, F, kind);
[~,s2] = qc_norm(C+0.99*F, F, kind);

zeta_bound = 2*max(abs(eig(F)))*max(abs(eig(C)))/min([abs(s1(2)), abs(s2(1))]);
fprintf('Absence bound for bang-bang: %s\n', zeta_bound);

%u0 = 0.5*ones(n,1);
%u0 = 1/m*ones(m,1);
dd = num2cell(readmatrix(sprintf('../njp_aroqa/dat/data_%s_norm%s_zeta%.0e.txt', parm, kind, zeta), 'Delimiter', ' '),1);
[u0,~,~,~] = deal(dd{:});

u = u0;
[u, ~, xe] = qc_optimize(C, F, u0, x0, kind, 'Zeta', zeta, 'Iterations', 300, 'Gamma', 2.2, 'Alpha', 0.9);

[J, hamiltonian, u_analytic, switching, lb, ub] = ...
    qc_analyse(u, C, F, x0, zeta, kind);

qc_save(parm, kind, J, lb, ub, zeta, u, u_analytic, switching, hamiltonian);
%writematrix([J, lb, ub, zeta], ...
%    sprintf('../dat/meta_%s_norm%s_zeta%.0e.txt', parm, kind, zeta), 'Delimiter', ' ');
%writematrix([u, u_analytic.', switching.', hamiltonian.'], ...
%    sprintf('../dat/data_%s_norm%s_zeta%.0e.txt', parm, kind, zeta), 'Delimiter', ' ');

subplot(2,1,1);
hold on
plot(linspace(0,1,size(u,1)), u);
plot(linspace(0,1,size(u,1)), u_analytic);
legend({'u', 'u_ana'});
ylim([0,1]);
hold off
subplot(2,1,2);
hold on
fill([0,1,1,0],[ub,ub,lb,lb],'k');
alpha(.25);
plot(linspace(0,1,size(u,1)), hamiltonian)
plot(linspace(0,1,size(u,1)), switching);
legend({'area', 'H', 'S'});
ylim([-3,4]);
hold off

