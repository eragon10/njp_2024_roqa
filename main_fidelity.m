
n = 200;
parm = 'x8params';

[B,C,F,x0] = eval(sprintf('%s()', parm));

%all = find_config(parm);
%all = ["x8params_norm2_zeta2e-01";
%       "x8params_normF_zeta1e-01";
%       "x8params_norm2sh_zeta3e-01"];

all = [ "x8params_norm2_zeta2e-01";
        "x8params_normF_zeta1e-01";
        "x8params_normnone_zeta0e+00";
        "x8params_normqaoa_zeta0e+00"];
%        "x8params_norm2sh_zeta3e-01"];
   
inputs = {};
kinds  = {};
zetas  = {};

names  = {};

for i = 1:size(all)
    [~,~,~,zetas{i},inputs{i}, ~,~,~,kinds{i}] = qc_load(all(i)); names{i} = sprintf('run_%s_%e', kinds{i}, zetas{i});
end

names = [{'*'}, names];
Jm = [];
Fm = [];


errors = {};
for i = 1:10
   errors{i} = min(1, max(-1, randn(10, 1)));
end

%errors{i+1} = (-1).^(1:10).';
%errors{i+2} = -1*(-1).^(1:4).';

for bb = 0.1:0.1:1.0
    er = {};
    for i = 1:size(bb,2)
        er{i} = bb*errors{i};
    end
    
    [Js, Fs, Bv] = qc_samples(C, F, x0, inputs, kinds, ...
       'Error', er, 'MaxStep', 3e-2);
   
    Jmeans = {};
    Fmins  = {};
    for i = 1:size(all)
        Jmeans{i} = mean(Js{i});
        Fmins{i}  = min(Fs{i});
    end

    Jm = [Jm; [sprintf('%f', bb), Jmeans]];
    Fm = [Fm; [sprintf('%f', bb), Fmins]];
end

writecell([names; Jm], sprintf('../njp_aroqa/dat/objmean_%s_full.csv', parm), 'Delimiter', ' ');
writecell([names; Fm], sprintf('../njp_aroqa/dat/fidmini_%s_full.csv', parm), 'Delimiter', ' ');

function [list] = find_config(paramname)
    list = [];
    fl = dir('../njp_aroqa/dat');
    for i = 1:size(fl)
        if and(contains(fl(i).name, paramname), contains(fl(i).name, 'data'))
           list = [list; string(sprintf('%s', fl(i).name(6:end-4)))]; 
        end
    end
end