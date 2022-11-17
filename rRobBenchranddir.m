function outstruct = rRobBenchranddir(args)

% r-Robustness Determination Benchmark
% Written by James Usevitch
%
% This function compares the performance of several algorithms on determining
% the maximal r-robustness of random digraphs.
%
% args.nvec : vector of n values to test. Column vector.
% args.pvec : vector of edge probability values to test. Row vector.
% args.test1robust : Set to 1 to ensure the graphs are at least 1 robust before testing
%                    them. Leave blank or set to 0 otherwise.
% args.repeatranddir : number of times to repeat the benchmark for each n and p
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

pool = gcp('nocreate')

% Create a pool if one isn't running already
if isempty(pool)
    %     pool = parpool(11); % or 12
    pool = parpool('local_Copy',args.poolsize);
end

addAttachedFiles(pool,{'../randDigraph.m'}); % Change to specify path to ErdRen.m file if necessary

% Create n vector (number of nodes in the graph)

nvec = args.nvec; % Make sure this is a column vector

% Erdos Renyi

% Vector of probabilities
% pvec = .3:.1:.8;
% pvec = pvec';

pvec = args.pvec;

algs = args.algs;

% Long vector of all n and p combinations; e.g. [n1 p1; n1 p2; ... n1 pend; n2 p1; ...]
% Makes it possible for the parfor loop to split each scenario between the
% workers more efficiently so no one core gets bogged down by calculating
% all graphs with a large n value or a large p value
nones = [ones(length(pvec),1) zeros(length(pvec),1)];
pones = [zeros(length(nvec),1) ones(length(nvec),1)];
totvec = kron(nvec,nones) + kron(pones,pvec); % Vector of all [n p] combos

if isfield(args,'repeatranddir') && args.repeatranddir > 1
    % Repeat the tests by duplicating the totvec as many times as specified
    % by args.repeat
    totvec = repmat(totvec,args.repeat,1);
    totvec = sortrows(totvec,[1 2]);
end


parfor ii=1:1:length(totvec)
    
    n = totvec(ii,1);
    p = totvec(ii,2);
    
    A = randDigraph(n,p);
    L = diag(A*ones(n,1)) - A;
    
    % Ensure robustness is at least 1, if args.test1robust == 1
    % THIS PART TO BE ADDED LATER USING THE TUTTE TREE TEST
    
    % LeBlanc algorithm
    if algs(1) == 1
        tic
        Lstruct = DetermineRobustness(struct('A',A,'smax',1));
        time1 = toc; % Remember to divide this by 2 since you check swapped S1 / S2 sets
        
        disp([newline 'Completed Random Directed Graph LeBlanc algorithm for n=' num2str(n) ' and p=' num2str(p) newline]); % char(10) makes a new line
        
        LB = [LB; [time1 Lstruct.r n p]];
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
        
        intlin = [intlin; [time2 istruct.maxr n p]];
        
        intlinwm = [intlinwm; [time3 iwmstruct.maxr n p]]
        
        disp(['Completed n=' num2str(n) ' and p=' num2str(p)])
    end
    
    testedMatrices = [testedMatrices; {L 'randdir' n p}];
end

outstruct.LB = LB;
outstruct.intlin = intlin;
outstruct.intlinwm = intlinwm;
outstruct.testedMatrices = testedMatrices;


disp(['Random Directed Graph analysis done' newline])



end