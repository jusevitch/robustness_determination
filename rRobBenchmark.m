% r-Robustness Determination Benchmark
% Written by James Usevitch
%
% WARNING: THIS SCRIPT WILL NOT RUN VERY FAST. Use the function
%   "rRobBenchmark2.m" for MUCH faster results, since functions are faster
%   than scripts in general.
%
% This script compares the performance of several algorithms on determining
% the maximal r-robustness of digraphs:
%   > (slightly modified) LeBlanc algorithm, with maximum value of s equal to
%       1
%   > General Integer Linear Program formulation by Usevitch and Panagou
%       (possibly nonempty S0 sets) using intlinprog
%   > General ILP using Gurobi
%   > General ILP using intlinprog using the two-set case as a warm start
%   > General ILP using Gurobi using the two-set case as a warm start
%
% NOTE: THIS SCRIPT MUST BE RUN ON MATLAB 2018a OR LATER. Any earlier
% version of Matlab does not allow you to specify an initial warm start
% point for intlinprog, and will throw an error.

clear all
clc

% Create structures to store metrics

finalLB = struct();
finalintlin = struct();
finalintlinwm = struct();
finalgurobi = struct();
finalgurobiwm = struct();

LBmetrics = [];
intlinmetrics = [];
intlinwm = [];

testGurobi = false; % Set to true to test with Gurobi
gurobimetrics = [];
gurobiwm = [];

testedMatrices = {}; % Initialize Matrix of tested matrices

% Part 1: LeBlanc vs ILP Algorithms --------------------------------------

% Open parpool - Test LeBlanc and intlinprog, save all matrices

pool = gcp('nocreate')

% Create a pool if one isn't running already
if isempty(pool)
    %     pool = parpool(11); % or 12
    pool = parpool('local_Copy',8   ); % For TV computer
end

addAttachedFiles(pool,{'../ErdRen.m'}); % Change to specify path to ErdRen.m file if necessary

% Create n vector (number of nodes in the graph)

nvec = 5:1:15;
nvec = nvec';



% Erdos Renyi

% Vector of probabilities
% pvec = .3:.1:.8;
% pvec = pvec';

pvec = [.3; .5; .8];

% Long vector of all n and p combinations; e.g. [n1 p1; n1 p2; ... n1 pend; n2 p1; ...]
% Makes it possible for the parfor loop to split each scenario between the
% workers more efficiently so no one core gets bogged down by calculating
% all graphs with a large n value or a large p value
nones = [ones(length(pvec),1) zeros(length(pvec),1)];
pones = [zeros(length(nvec),1) ones(length(nvec),1)];
totvec = kron(nvec,nones) + kron(pones,pvec); % Vector of all [n p] combos


parfor ii=1:1:length(totvec)
    
    n = totvec(ii,1);
    p = totvec(ii,2);
    
    A = ErdRen(n,p);
    L = diag(A*ones(n,1)) - A;
    
    % LeBlanc algorithm
    tic
    Lstruct = DetermineRobustness(struct('A',A,'smax',1));
    time1 = toc; % Remember to divide this by 2 since you check swapped S1 / S2 sets
    
    disp([newline 'Completed Erdos-Renyi LeBlanc algorithm for n=' num2str(n) ' and p=' num2str(p) newline]); % char(10) makes a new line
    
    % intlinprog algorithm
    tic
    istruct = rRobustGnl(struct('L',L));
    time2 = toc;
    
    % intlinprog warm start
    tic
    iwmstruct = rRobustWm(struct('L',L));
    time3 = toc;
    
    testedMatrices = [testedMatrices; {L 'Erdos' n p}];
    
    LBmetrics = [LBmetrics; [time1 Lstruct.r n p]];
    
    intlinmetrics = [intlinmetrics; [time2 istruct.maxr n p]];
    
    intlinwm = [intlinwm; [time3 iwmstruct.maxr n p]]
    
    disp(['Completed n=' num2str(n) ' and p=' num2str(p)])
end

finalLB.Erdos = LBmetrics;
finalintlin.Erdos = intlinmetrics;
finalintlinwm.Erdos = intlinwm;

disp(['Erdos-Renyi analysis done' newline])




% Complete graphs

LBmetrics = [];
intlinmetrics = [];
intlinwm = [];

time1 = 0;
time2 = 0;
time3 = 0;

% Types of analysis: 1 - Leblanc, 2 - intlinprog, 3 - intlinprogwarm

nones = [ones(3,1) zeros(3,1)];
aones = [zeros(length(nvec),1) ones(length(nvec),1)];
totvec = kron(nvec,nones) + kron(aones,[1;2;3]);

parfor ii=1:1:length(totvec)
    
    n = totvec(ii,1);
    
    L = makegraph(struct('n',n,'type','complete'));
    A = diag(diag(L)) - L;
    
    testedMatrices = [testedMatrices; {L 'complete' n 0}];
    
    
    if totvec(ii,2) == 1
        % LeBlanc algorithm
        tic
        Lstruct = DetermineRobustness(struct('A',A,'smax',1));
        time1 = toc; % No longer need to divide by 2
        LBmetrics = [LBmetrics; [time1 Lstruct.r n]];
        disp([newline 'Completed LeBlanc algorithm for complete graphs with n=' num2str(n) newline]);
        
    elseif totvec(ii,2) == 2
        % intlinprog algorithm
        tic
        istruct = rRobustGnl(struct('L',L));
        time2 = toc;
        intlinmetrics = [intlinmetrics; [time2 istruct.maxr n]];
        
    elseif totvec(ii,2) == 3
        % intlinprog warm start
        tic
        iwmstruct = rRobustWm(struct('L',L));
        time3 = toc;
        intlinwm = [intlinwm; [time3 iwmstruct.maxr n]]
        
    else
        error('Something bad happened in the complete graph section...');
    end
    
    disp(['Completed n=' num2str(n)])
end

finalLB.complete = LBmetrics;
finalintlin.complete = intlinmetrics;
finalintlinwm.complete = intlinwm;

disp(['Complete graph analysis done' newline])



% k-Circulant Undirected

LBmetrics = [];
intlinmetrics = [];
intlinwm = [];

time1 = 0;
time2 = 0;
time3 = 0;

maxk = 0;

% Define kvec out here if you want it to be the same regardless of n
% kvec = [3 5 10];

parfor ii=1:1:length(nvec)
    
    n = nvec(ii);
    
    % Use this section if you want to auto-generate all possible values of
    % k
    % If n is even, make sure that max k is (n/2) - 1
    if mod(n,2) == 0
        maxk = n/2 -1;
    else
        maxk = floor(n/2);
    end
    kvec = 2:3:maxk;
    %     kvec = 2:2:maxk; % Perhaps use this if n gets too big
    for k=1:1:length(kvec)
        
        L = makegraph(struct('n',n,'k',kvec(k),'type','kundir'));
        A = diag(diag(L)) - L;
        
        % LeBlanc algorithm
        tic
        Lstruct = DetermineRobustness(struct('A',A,'smax',1));
        time1 = toc; % Remember to divide this by 2 since you check swapped S1 / S2 sets
        
        disp([newline 'Completed k-circulant Undirected LeBlanc algorithm for n=' num2str(n) ' and k=' num2str(kvec(k)) newline]);
        
        % intlinprog algorithm
        tic
        istruct = rRobustGnl(struct('L',L));
        time2 = toc;
        
        % intlinprog warm start
        tic
        iwmstruct = rRobustWm(struct('L',L));
        time3 = toc;
        
        testedMatrices = [testedMatrices; {L 'kundir' n kvec(k)}];
        
        LBmetrics = [LBmetrics; [time1 Lstruct.r n kvec(k)]];
        
        intlinmetrics = [intlinmetrics; [time2 istruct.maxr n kvec(k)]];
        
        intlinwm = [intlinwm; [time3 iwmstruct.maxr n kvec(k)]]
        
        
    end
    disp(['Completed n=' num2str(n) ' and all values of k'])
end

finalLB.kundir = LBmetrics;
finalintlin.kundir = intlinmetrics;
finalintlinwm.kundir = intlinwm;

disp(['k-Circulant undirected graphs done' newline])




% k-Circulant Directed

LBmetrics = [];
intlinmetrics = [];
intlinwm = [];

time1 = 0;
time2 = 0;
time3 = 0;

parfor ii=1:1:length(nvec)
    
    n = nvec(ii);
    
    maxk = n-2; % Technically can be n-1, but k = n-1 is a complete graph which was already tested
    kvec = 2:3:maxk;
    %     kvec = 2:2:maxk; % Perhaps use this if n gets too big
    for k=1:1:length(kvec)
        
        L = makegraph(struct('n',n,'k',kvec(k),'type','kdir'));
        A = diag(diag(L)) - L;
        
        % LeBlanc algorithm
        tic
        Lstruct = DetermineRobustness(struct('A',A,'smax',1));
        time1 = toc; % Remember to divide this by 2 since you check swapped S1 / S2 sets
        disp([newline 'Completed k-circulant Directed LeBlanc algorithm for n=' num2str(n) ' and k=' num2str(kvec(k)) newline]);
        
        % intlinprog algorithm
        tic
        istruct = rRobustGnl(struct('L',L));
        time2 = toc;
        
        % intlinprog warm start
        tic
        iwmstruct = rRobustWm(struct('L',L));
        time3 = toc;
        
        testedMatrices = [testedMatrices; {L 'kdir' n kvec(k)}];
        
        LBmetrics = [LBmetrics; [time1 Lstruct.r n kvec(k)]];
        
        intlinmetrics = [intlinmetrics; [time2 istruct.maxr n kvec(k)]];
        
        intlinwm = [intlinwm; [time3 iwmstruct.maxr n kvec(k)]]
        
        
    end
    
    disp(['Completed n=' num2str(n) ' for all k values'])
end

finalLB.kdir = LBmetrics;
finalintlin.kdir = intlinmetrics;
finalintlinwm.kdir = intlinwm;

disp(['k-Circulant directed graphs done.' newline])

if testGurobi
    disp('Moving to Gurobi implementations...')
    
    % Close parpool, test Gurobi implementations on all saved matrices
    
    delete(pool)
    
    % Gurobi implementation
    
    time1 = 0;
    time2 = 0;
    
    finalgurobi.Erdos = [];
    finalgurobi.complete = [];
    finalgurobi.kundir = [];
    finalgurobi.kdir = [];
    
    finalgurobiwm.Erdos = [];
    finalgurobiwm.complete = [];
    finalgurobiwm.kundir = [];
    finalgurobiwm.kdir = [];
    
    totmat = size(testedMatrices,1);
    disp([newline 'Total matrices to test:' ' ' num2str(totmat) newline])
    
    for ii=1:1:size(testedMatrices,1)
        
        
        cell = testedMatrices(ii,:);
        
        L = cell{1};
        n = cell{3};
        
        if strcmp(cell{2},'Erdos')
            p = cell{4};
            
            tic
            gstruct = rRobustGnl_gurobi(struct('L',L));
            time1 = toc;
            
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobi.Erdos = [finalgurobi.Erdos; [time1 gstruct.maxr n p]];
            
            finalgurobiwm.Erdos = [finalgurobiwm.Erdos; [time2 gstructwm.maxr n p]];
            
            
        elseif strcmp(cell{2},'complete')
            
            tic
            gstruct = rRobustGnl_gurobi(struct('L',L));
            time1 = toc;
            
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobi.complete = [finalgurobi.complete; [time1 gstruct.maxr n]];
            
            finalgurobiwm.complete = [finalgurobiwm.complete; [time2 gstructwm.maxr n]];
            
            
            
        elseif strcmp(cell{2},'kundir')
            k = cell{4};
            
            tic
            gstruct = rRobustGnl_gurobi(struct('L',L));
            time1 = toc;
            
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobi.kundir = [finalgurobi.kundir; [time1 gstruct.maxr n k]];
            
            finalgurobiwm.kundir = [finalgurobiwm.kundir; [time2 gstructwm.maxr n k]];
            
            
        elseif strcmp(cell{2},'kdir')
            k = cell{4};
            
            tic
            gstruct = rRobustGnl_gurobi(struct('L',L));
            time1 = toc;
            
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobi.kdir = [finalgurobi.kdir; [time1 gstruct.maxr n k]];
            
            finalgurobiwm.kdir = [finalgurobiwm.kdir; [time2 gstructwm.maxr n k]];
            
        end
        
        
        if mod(ii,10) == 0
            disp([num2str(ii) ' / ' num2str(totmat) ' matrices analyzed'])
        end
        
    end
    
    disp('Gurobi analysis done')
end

disp('Part 1 done. Moving on to Part 2...')



% Part 2: ILP Algorithms on Matrices of Known Robustness ------------------

% Do NOT save matrices


% nvec = [30; 40; 50; 75; 100];

% Complete Graphs





% k-Circulant Undirected Graphs




% k-Circulant directed Graphs



% k-Platoon




% Test Example
% args = struct();
%
% args.n = 20;
% args.k = 5;
% args.p = .3;
% args.type = 'complete';

% L = makegraph(args);
%
% args2 = struct('L',L);
%
% results1 = rRobustGnl(args2)
% results2 = rRobustWm(args2)
% results3 = rRobustGnl_gurobi(args2)
% results4 = rRobustWm_gurobi(args2)





