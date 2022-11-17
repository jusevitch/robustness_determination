function robustpaperdata(folder_name)

data_path = ['~/james/PhDpriv/Graph Theory/robust_determ/data/' folder_name];

if exist(data_path,'dir') ~= 7
    mkdir(data_path);
end

% rng('shuffle')

% DetermineRobustness portion

args1.nvec = [7:15]';
% args1.nvec = kron(args1.nvec, ones(100,1));
args1.nvec = kron(args1.nvec, ones(100,1));
args1.pvec = [.3;.5;.8];
args1.kvec = [3;4;5];
% args1.kvec = [3;4;5];
args1.MaxTime = 2*10^3;

args1.poolsize = 15;
args1.graphs = [1;1;1;1];
%   Algorithm groups are [LB, DetRob, rRobustGnl, rRobustGnl2, rsGnl2, Fmax, rRoblower, rRobupper]
args1.algs = [1;1;0;1;1;0;1;1];
% args1.MaxTime = 10;

results_7_15 = rsBench4(args1);
save([data_path '/results_7_15']);


% No DetermineRobustness portion

args2 = args1;
args2.nvec = [17:2:25]';
args2.nvec = kron(args2.nvec, ones(100,1));
args2.pvec = [.3; .5; .8];
args2.kvec = [3;4;5];

args2.algs(1:2) = [0;0]; % Don't include the DetermineRobustness algorithms. They'll take forever.
args2.algs(6) = 0; % Excludes maxF.m

results_17_25 = rsBench4(args2);

save([data_path '/results_17_25']);



% No Fmax portion

% args3 = args2;
% args3.nvec = [35:36]';
% args3.pvec = [.3; .5; .8];
% args3.kvec = [10];
% args3.print = 1;
% % args3.algs(4:5) = [0;0];
% args3.algs = [0;0;0;1;0;0;0;0];
% args3.MaxTime = 5;
% 
% results_35_50 = rsBench3(args3);
% 
% save([data_path '/results_35_50']);

end