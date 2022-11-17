function outstruct = rRobBenchkoutrand(args)

% r-Robustness Determination Benchmark
% Written by James Usevitch
%
% This function compares the performance of several algorithms on determining
% the maximal r-robustness of k-out random digraphs.
% See Bollobas 2001, "Random Graphs", Chapter 2.
%
% args.nvec : vector of n values to test. Column vector.
% args.kvec : vector of k values to test. Row vector.
% args.test1robust : Set to 1 to ensure the graphs are at least 1 robust before testing
%                    them. Leave blank or set to 0 otherwise.
% args.repeatkoutrand : number of times to repeat the benchmark for each n and p
%               pair. Integer value greater than or equal to 1.
% args.poolsize : number of workers to include in the parpool
% args.algs : Set entry to 1 to test the algorithm group, 0 otherwise. Algorithm groups are [LB and intlinprog and intlinprogwm; Gurobi and Gurobiwm]
%
% NOTE: THIS FUNCTION MUST BE RUN ON MATLAB 2018a OR LATER. Any earlier
% version of Matlab does not allow you to specify an initial warm start
% point for intlinprog, and will throw an error.

LB = [];
intlin = [];
intlinwm = [];

testedMatrices = {}; % Initialize Matrix of tested matrices

% Part 1: LeBlanc vs ILP Algorithms --------------------------------------

% Open parpool - Test LeBlanc and intlinprog, save all matrices

pool = gcp('nocreate');

% Create a pool if one isn't running already
if isempty(pool)
    %     pool = parpool(11); % or 12
    pool = parpool('local_Copy',args.poolsize);
end

addAttachedFiles(pool,{'../koutRandDigraph.m'}); % Change to specify path to ErdRen.m file if necessary

% Create n vector (number of nodes in the graph)

nvec = args.nvec; % Make sure this is a column vector

kvec = args.kvec;

algs = args.algs;

% Long vector of all n and k combinations; e.g. [n1 k1; n1 k2; ... n1 kend; n2 k1; ...]
% Makes it possible for the parfor loop to split each scenario between the
% workers more efficiently so no one core gets bogged down by calculating
% all graphs with a large n value or a large k value
nones = [ones(length(kvec),1) zeros(length(kvec),1)];
pones = [zeros(length(nvec),1) ones(length(nvec),1)];
totvec = kron(nvec,nones) + kron(pones,kvec); % Vector of all [n p] combos

if isfield(args,'repeatkoutrand') && args.repeatkoutrand > 1
    % Repeat the tests by duplicating the totvec as many times as specified
    % by args.repeat
    totvec = repmat(totvec,args.repeat,1);
    totvec = sortrows(totvec,[1 2]);
end


parfor ii=1:1:length(totvec)
    
    n = totvec(ii,1);
    k = totvec(ii,2);
    
    A = koutRandDigraph(n,k);
    L = diag(A*ones(n,1)) - A;
    
    % Ensure robustness is at least 1, if args.test1robust == 1
    % THIS PART TO BE ADDED LATER USING THE TUTTE TREE TEST
    
    % LeBlanc algorithm
    if algs(1) == 1
        tic
        Lstruct = DetermineRobustness(struct('A',A,'smax',1));
        time1 = toc; % Remember to divide this by 2 since you check swapped S1 / S2 sets
        
        disp([newline 'Completed k-out LeBlanc algorithm for n=' num2str(n) ' and k=' num2str(k) newline]); % char(10) makes a new line
        
        LB = [LB; [time1 Lstruct.r n k]];
    end
    
    % intlinprog algorithm
    if algs(2) == 1
        tic
        istruct = rRobustGnl(struct('L',L));
        time2 = toc;
        
        % intlinprog warm start
        tic
        iwmstruct = rRobustWm(struct('L',L));
        time3 = toc;
        
        intlin = [intlin; [time2 istruct.maxr n k]];
        
        intlinwm = [intlinwm; [time3 iwmstruct.maxr n k]]
        
        disp(['Completed n=' num2str(n) ' and k=' num2str(k)])
    end
    
    testedMatrices = [testedMatrices; {L 'koutrand' n k}];
end

outstruct.LB = LB;
outstruct.intlin = intlin;
outstruct.intlinwm = intlinwm;
outstruct.testedMatrices = testedMatrices;


disp(['k-out RD Graphs analysis done' newline])



end