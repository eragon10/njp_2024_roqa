
n = 200;
parm = 'x8params';

[B,C,F,x0] = eval(sprintf('%s()', parm));

all = find_config(parm);

inputs = {};
kinds  = {};
zetas  = {};

for i = 1:size(all)
    [~,~,~,zetas{i},inputs{i}, ~,~,~,kinds{i}] = qc_load(all(i));
end

er = {};
bb = 0.5;

for i = 1:100
   er{i} = min(bb, max(-bb, randn(20, 1)));
end
%e = readmatrix('error.txt');
%e = readmatrix('x3err.txt');
%e = min(bb, max(-bb, sum([inputs{:}],2)/size(inputs,2)));
%e = min(bb, max(-bb, randn(20, 1)));

% C, F, x0, us, opts, bound, samples
[Js, Fs, Bv] = qc_samples(C, F, x0, inputs, kinds, ...
   'Error', er, 'MaxStep', 3e-2); % 0.2*(-1).^(1:10).'

for i = 1:size(all)
    fprintf(' => |norm:%s\t|zeta:%f| variance:%f | fidelity:%f \n', kinds{i}, zetas{i}, mean(Js{i}), min(Fs{i}));
end


function [list] = find_config(paramname)
    list = [];
    fl = dir('../dat');
    for i = 1:size(fl)
        if and(contains(fl(i).name, paramname), contains(fl(i).name, 'data'))
           list = [list; string(sprintf('%s', fl(i).name(6:end-4)))]; 
        end
    end
end