function outstruct = rsErdos(args)

% (r,s)-Robustness Determination Benchmark
% Written by James Usevitch
%
% This function compares the performance of several algorithms on determining
% the max values of r and s for which graphs are (r,s)-robust, and the max
% values of F for which graphs are (F+1,F+1)-robust
%
% args.nvec : vector of n values to test. Column vector.
% args.pvec : vector of edge probability values to test. Row vector.
% args.test1robust : Set to 1 to ensure the graphs are at least 1 robust before testing
%                    them. Leave blank or set to 0 otherwise.
% args.repeatErdos : number of times to repeat the benchmark for each n and p
%               pair. Integer value greater than or equal to 1.
% args.poolsize : number of workers to include in the parpool
% args.algs : Set entry to 1 to test the algorithm group, 0 otherwise. 
%   Algorithm groups are [LB, DetRob, rRobustGnl, rRobustGnl2, rsGnl2, Fmax] (intlinprog refers to max (r,s) robustness) 


% Determine algorithms to test
if ~isfield(args,'algs')
    algs = [1;1;1;1;1;1]; % Test all algorithms
else
    algs = args.algs;
end

if isempty(find(algs))
    error('args.algs : No algorithm specified')
end

% Create structures to store metrics

LB = [];
intlin = [];
maxf = [];


testedMatrices = {}; % Initialize Matrix of tested matrices

% Part 1: LeBlanc vs ILP Algorithms --------------------------------------

% Open parpool - Test LeBlanc and intlinprog, save all matrices

pool = gcp('nocreate')

% Create a pool if one isn't running already
if isempty(pool)
    %     pool = parpool(11); % or 12
    pool = parpool('local_Copy',args.poolsize);
end

addAttachedFiles(pool,{'../ErdRen.m'}); % Change to specify path to ErdRen.m file if necessary

% Create n vector (number of nodes in the graph)

nvec = args.nvec; % Make sure this is a column vector

% Erdos Renyi

pvec = args.pvec;

% Long vector of all n and p combinations; e.g. [n1 p1; n1 p2; ... n1 pend; n2 p1; ...]
% Makes it possible for the parfor loop to split each scenario between the
% workers more efficiently so no one core gets bogged down by calculating
% all graphs with a large n value or a large p value
nones = [ones(length(pvec),1) zeros(length(pvec),1)];
pones = [zeros(length(nvec),1) ones(length(nvec),1)];
totvec = kron(nvec,nones) + kron(pones,pvec); % Vector of all [n p] combos

if isfield(args,'repeatErdos') && args.repeatErdos > 1
    % Repeat the tests by duplicating the totvec as many times as specified
    % by args.repeat
    totvec = repmat(totvec,args.repeat,1);
    totvec = sortrows(totvec,[1 2]);
end

% Create master list of total matrices

parfor ii=1:1:length(totvec)
    n = totvec(ii,1);
    p = totvec(ii,2);
    
    A = ErdRen(n,p);
    L = diag(A*ones(n,1)) - A;
    
    % Ensure robustness is at least 1, if args.test1robust == 1
    if isfield(args,'test1robust') && args.test1robust
        eigvals = sort(eig(L));
        if eigvals(2) <= 1e-11
            iters = 1;
            while eigvals(2) <= 1e-11 % Effectively zero
                disp(['Attempting to find nonsingular A. Try ', ' ', num2str(iters)])
                A = ErdRen(n,p);
                L = diag(A*ones(n,1)) - A;
                eigvals = sort(eig(L));
                iters = iters + 1;
            end
        end
    end
    
    testedMatrices = [testedMatrices; {L 'Erdos' n p A}];
    
end



for ii=1:1:size(testedMatrices,1)
    
    L = testedMatrices{ii,1};
    A = testedMatrices{ii,end};
    n = testedMatrices{ii,3};
    p = testedMatrices{ii,4};
    
    % LB algorithm -- smax == 1
    if algs(1) == 1
        Lstruct = DetermineRobustness(struct('A',A));
        time1 = Lstruct.time;
        
        
        disp([newline 'Completed Erdos-Renyi LeBlanc algorithm for n=' num2str(n) ' and p=' num2str(p) newline]); % char(10) makes a new line
        
        LB = [LB; [time1 Lstruct.r n p Lstruct.s]];
    end
    
    if algs(2) == 1
        
        
        
    end
    
    % Max (r,s)-robustness
    if algs(3) == 1
        istruct = rsGnl2(struct('L',L));
        time2 = istruct.totaltime;
        
        intlin = [intlin; [time2 istruct.r n p istruct.s]];
        
        disp(['Completed (r,s) intlin for n=' num2str(n) ' and p=' num2str(p)])
    end
    
    if algs(4) == 1
        fstruct = maxF(struct('L',L));
        time3 = fstruct.time;
        
        maxf = [maxf; [time3 fstruct.Fhat n p]];
        disp(['Completed Fmax for n=' num2str(n) ' and p=' num2str(p)])
    end
end

outstruct.LB = LB;
outstruct.intlin = intlin;
outstruct.intlinwm = intlinwm;
outstruct.testedMatrices = testedMatrices;


disp(['Erdos-Renyi analysis done' newline])





end