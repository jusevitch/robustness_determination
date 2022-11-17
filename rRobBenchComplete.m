function outstruct = rRobBenchComplete(args)

% r-Robustness Determination Benchmark - Complete
% Written by James Usevitch
%
% This function compares the performance of several algorithms on determining
% the maximal r-robustness of complete graphs
%
% args.nvec : vector of n values to test. Column vector.
% args.poolsize : number of workers to include in the parpool
%
% NOTE: THIS FUNCTION MUST BE RUN ON MATLAB 2018a OR LATER. Any earlier
% version of Matlab does not allow you to specify an initial warm start
% point for intlinprog, and will throw an error.


% Create structures to store metrics

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
    pool = parpool('local_Copy', args.poolsize); % For TV computer
end

% Create n vector (number of nodes in the graph)

nvec = args.nvec; % Make sure this is a column vector


% Complete graphs

LB = [];
intlin = [];
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
        LB = [LB; [time1 Lstruct.r n]];
        disp([newline 'Completed LeBlanc algorithm for complete graphs with n=' num2str(n) newline]);
        
    elseif totvec(ii,2) == 2
        % intlinprog algorithm
        tic
        istruct = rRobustGnl(struct('L',L));
        time2 = toc;
        intlin = [intlin; [time2 istruct.maxr n]];
        
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

% Every complete matrix is repeated 3 times in testedMatrices, so reduce
% this to only once per matrix

idxtemp = 1:3:size(testedMatrices,1);
idxtemp = idxtemp';
testedMatrices = testedMatrices(idxtemp,:);

outstruct.LB = LB;
outstruct.intlin = intlin;
outstruct.intlinwm = intlinwm;
outstruct.testedMatrices = testedMatrices;

disp(['Complete graph analysis done' newline])

end