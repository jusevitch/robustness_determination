function outcell = gurobiGraphs(args)

% Make a bunch of matrices to test for the Gurobi algorithms
% args.nvec :
% args.algs : Set entry to 1 to test the algorithm group, 0 otherwise. Algorithm groups are [LB; intlinprog and intlinprogwm; Gurobi and Gurobiwm]
% args.graphs : Set entry to 1 to test the graph type, 0 otherwise. Graph types are (see makegraph.m function for details):
%               [Erdos, complete, kundir, kdir, randdir, kinrand, koutrand]


% Open parpool - Test LeBlanc and intlinprog, save all matrices


outcell = {};

% Parfor, if it makes things go faster. I don't think it does though.
% pool = gcp('nocreate')
%
% % Create a pool if one isn't running already
% if isempty(pool)
%     %     pool = parpool(11); % or 12
%     pool = parpool('local_Copy',args.poolsize);
% end
%
% addAttachedFiles(pool,{'../ErdRen.m'}); % Change to specify path to ErdRen.m file if necessary


% Erdos graphs

if args.graphs(1)
    nvec = args.nvec; % Make sure this is a column vector
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
    
    testedMatrices = {}; % Initialize Matrix of tested matrices
    
    for ii=1:1:length(totvec)
        %     parfor ii=1:1:length(totvec)
        
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
        
        testedMatrices = [testedMatrices; {L 'Erdos' n p}];
        
    end
    
    outcell = [outcell; testedMatrices];
    
end


% Complete graphs

if args.graphs(2)
    testedMatrices = {};
    
    nvec = args.nvec;
    
    for ii=1:1:length(nvec)
        
        n = nvec(ii);
        L = makegraph(struct('n',n,'type','complete'));
        A = diag(diag(L)) - L;
        
        testedMatrices = [testedMatrices; {L 'complete' n 0}];
    end
    
    outcell = [outcell; testedMatrices];
    
end


% kundir graphs

if args.graphs(3)
    
    % TBD
    
end


% kdir graphs

if args.graphs(4)
    
    % TBD
    
end


% randdir graphs

if args.graphs(5)
    
    testedMatrices = {};
    nvec = args.nvec;
    pvec = args.pvec;
    
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
    
    for ii=1:1:length(totvec)
        n = totvec(ii,1);
        p = totvec(ii,2);
        
        A = randDigraph(n,p);
        L = diag(A*ones(n,1)) - A;
        
        testedMatrices = [testedMatrices; {L 'randdir' n p}];
    end
    
    outcell = [outcell; testedMatrices];
end


% kinrand graphs

if args.graphs(6)
    
    testedMatrices = {};
    
    nvec = args.nvec; % Make sure this is a column vector
    
    kvec = args.kvec;
    
    nones = [ones(length(kvec),1) zeros(length(kvec),1)];
    pones = [zeros(length(nvec),1) ones(length(nvec),1)];
    totvec = kron(nvec,nones) + kron(pones,kvec); % Vector of all [n p] combos
    
    if isfield(args,'repeatkinrand') && args.repeatkinrand > 1
        % Repeat the tests by duplicating the totvec as many times as specified
        % by args.repeat
        totvec = repmat(totvec,args.repeat,1);
        totvec = sortrows(totvec,[1 2]);
    end
    
    for ii=1:1:length(totvec)
        n = totvec(ii,1);
        k = totvec(ii,2);
        
        A = kinRandDigraph(n,k);
        L = diag(A*ones(n,1)) - A;
        
        testedMatrices = [testedMatrices; {L 'kinrand' n k}];
        
        
    end
    
    outcell = [outcell; testedMatrices];
end



% koutrand graphs

if args.graphs(7)
    testedMatrices = {};
    
    nvec = args.nvec; % Make sure this is a column vector
    
    kvec = args.kvec;
    
    nones = [ones(length(kvec),1) zeros(length(kvec),1)];
    pones = [zeros(length(nvec),1) ones(length(nvec),1)];
    totvec = kron(nvec,nones) + kron(pones,kvec); % Vector of all [n p] combos
    
    if isfield(args,'repeatkoutrand') && args.repeatkoutrand > 1
        % Repeat the tests by duplicating the totvec as many times as specified
        % by args.repeat
        totvec = repmat(totvec,args.repeat,1);
        totvec = sortrows(totvec,[1 2]);
    end
    
    for ii=1:1:length(totvec)
        
        n = totvec(ii,1);
        k = totvec(ii,2);
        
        A = koutRandDigraph(n,k);
        L = diag(A*ones(n,1)) - A;
        
        testedMatrices = [testedMatrices; {L 'koutrand' n k}];
    end
    
    outcell = [outcell; testedMatrices];
end

end