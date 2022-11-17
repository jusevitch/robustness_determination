function outstruct = rRobBenchkdir(args)

% r-Robustness Determination Benchmark - k-Circulant Directed Graphs
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
% kvec = args.kvec;

% k-Circulant Undirected

LB = [];
intlin = [];
intlinwm = [];

time1 = 0;
time2 = 0;
time3 = 0;

maxk = 0;

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
        
        LB = [LB; [time1 Lstruct.r n kvec(k)]];
        
        intlin = [intlin; [time2 istruct.maxr n kvec(k)]];
        
        intlinwm = [intlinwm; [time3 iwmstruct.maxr n kvec(k)]]
        
        
    end
    
    disp(['Completed n=' num2str(n) ' for all k values'])
end

outstruct.LB = LB;
outstruct.intlin = intlin;
outstruct.intlinwm = intlinwm;
outstruct.testedMatrices = testedMatrices;

disp(['k-Circulant directed graphs done.' newline])


end