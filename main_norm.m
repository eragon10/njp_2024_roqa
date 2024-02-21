

kind = '2';

n = 1000;
o = 2;

q = zeros(n+2*o,1);

v = zeros(n+2*o,1);
su = zeros(n+2*o,1);
sl = zeros(n+2*o,1);


v1 = zeros(n+2*o,1);
s1 = zeros(n+2*o,1);


kind = '2';

parm = 'x8params';

[B,C,F,x0] = eval(sprintf('%s()', parm));


for i = 1:(n+2*o)
    q(i) = (i-o) / n;
    [v(i), subs] = qc_norm(C + q(i) * F, F, kind);
    sl(i) = subs(1);
    su(i) = subs(2); %0.5*sum(subs);
    
    %[v1(i), subs] = qc_norm(C + q(i) * F, F, '2sh_fd');
    %s1(i) = subs(2); %0.5*sum(subs);
    
end

writematrix([q, v, sl, su], sprintf('../dat/plot_%s_norm%s.txt', parm, kind), 'Delimiter', ' ');


plot(q, v, q, sl, q, su) %, q, v1, q, s1);
legend({'sigma', 'subdif_l', 'subdif_u' 'v1', 's1'});