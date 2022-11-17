function outstruct = rsBench3(args)

% args.nvec : vector of n values to test. Column vector.
% args.pvec : vector of edge probability values to test. Column vector.
% args.kvec : vector of k values to test. Column vector.
% args.test1robust : Set to 1 to ensure the graphs are at least 1 robust before testing
%                    them. Leave blank or set to 0 otherwise.
% args.poolsize : number of workers to include in the parpool
% args.graphs : Set entry to 1 to test the graph type, 0 otherwise. Graph types are (see makegraph.m function for details):
%               [Erdos, randdir, kinrand, koutrand]
%               *** p graphs go first, k graphs go last ***
% args.algs : Set entry to 1 to test the algorithm group, 0 otherwise.
%   Algorithm groups are [LB, DetRob, rRobustGnl, rRobustGnl2, rsGnl2, Fmax, rRoblower, rRobupper] (intlinprog refers to max (r,s) robustness)

%
% Check args conditions, create variables, open parpool
%

% Determine algorithms to test
if ~isfield(args,'algs')
    algs = [1;1;1;1;1;1;1;1]; % Test all algorithms
else
    algs = args.algs;
end

if ~isfield(args,'graphs')
    graphs = [1;1;1;1];
else
    graphs = args.graphs;
end

if isempty(find(algs))
    error('args.algs : No algorithm specified')
end

if ~isfield(args,'print')
    print = 0;
else
    print = args.print;
end

% Create structures to store metrics

% Erdos = [];
% randdir = [];
% kin = [];
% kout = [];

LB = [];
DetRob = [];
rRG = [];
rRG2 = [];
rs2 = [];
Fmax = [];
rlow = [];
rup = [];

testedMatrices = {}; % Initialize Matrix of tested matrices

% List of graph types
graphnames = {'Erdos', 'randdir', 'kinrand', 'koutrand'};
% algnames = {LB, DetRob, rRobustGnl, rRobustGnl2, rsGnl2, Fmax}

% Open parpool - Test LeBlanc and intlinprog, save all matrices

pool = gcp('nocreate')

% Create a pool if one isn't running already
if isempty(pool)
    %     pool = parpool(11); % or 12
    pool = parpool('local_Copy',args.poolsize);
end

% addAttachedFiles(pool,{'../ErdRen.m'}); % Change to specify path to ErdRen.m file if necessary


%
% Create ALL graphs
%

nvec = args.nvec; % Make sure this is a column vector
pvec = args.pvec;
kvec = args.kvec;
avec = find(args.algs); % Number representation of algorithms


% % Long vector of all n, p/k, and algorithm combinations.
% % This allows all situations to be split between the threads, preventing
% % long run times on one thread for a high value of n
% % n* = |{n1,...nend}|
% % pk* = |{pk1,...,pkend}|
% % a* = |{alg1,...,algend}|
% % For all combinations of n, pk, and alg:
% %   1st column n, 2nd column pk, last column alg #
% %   Repeat alg list (n*)(pk*) times.
% %   Form sublist pksub by repeating each entry of pk a* times. Repeat pksub
% %   n* times.
% %   Repeat each entry of n list (a*)(pk*) times.
% totvecp = zeros(length(nvec)*length(pvec)*length(avec),3);
% totveck = zeros(length(nvec)*length(kvec)*length(avec),3);
% 
% totvecp(:,3) = repmat(avec, length(nvec)*length(pvec),1);
% psub = repelem(pvec,length(avec));
% totvecp(:,2) = repmat(psub, length(nvec),1);
% totvecp(:,1) = repelem(nvec,length(pvec)*length(avec));
% 
% totveck(:,3) = repmat(avec, length(nvec)*length(kvec),1);
% ksub = repelem(kvec,length(avec));
% totveck(:,2) = repmat(ksub, length(nvec),1);
% totveck(:,1) = repelem(nvec,length(kvec)*length(avec));


% Long vector of all n and p combinations; e.g. [n1 p1; n1 p2; ... n1 pend; n2 p1; ...]
% Makes it possible for the parfor loop to split each scenario between the
% workers more efficiently so no one core gets bogged down by calculating
% all graphs with a large n value or a large p value
nones = [ones(length(pvec),1) zeros(length(pvec),1)];
pones = [zeros(length(nvec),1) ones(length(nvec),1)];
totvecp = kron(nvec,nones) + kron(pones,pvec); % Vector of all [n p] combos

% Long vector of all n and k combinations; e.g. [n1 k1; n1 k2; ... n1 kend; n2 k1; ...]
nones = [ones(length(kvec),1) zeros(length(kvec),1)];
kones = [zeros(length(nvec),1) ones(length(nvec),1)];
totveck = kron(nvec,nones) + kron(kones,kvec); % Vector of all [n p] combos

% Create master list of total matrices

lastpidx = 2; % index of the last graph requiring p parameter in graphs vector

for ii=1:1:length(graphs)
    graphtype = graphnames(ii);
    if ii <= lastpidx
        totvec = totvecp;
    else
        totvec = totveck;
    end
    
    parfor jj=1:1:length(totvec)
        n = totvec(jj,1);
        pk = totvec(jj,2); % Represents either p or k, depending on which graph is being treated
%         algnum = totvec(jj,3);
        
        if ii <= lastpidx
            L = makegraph(struct('n',n,'p',pk,'type',graphtype));
        else
            L = makegraph(struct('n',n,'k',pk,'type',graphtype));
        end
        
        % Ensure robustness is at least 1, if args.test1robust == 1
        % To be added later
        
        testedMatrices = [testedMatrices; {L ii n pk}]; % The second index is the graphtype number
    end
end

testedMatrices = [repelem(testedMatrices,length(avec),1) num2cell(repmat(avec,length(testedMatrices),1))];



%
% Perform Algorithm testing
% Algorithm groups are [LB, DetRob, rRobustGnl, rRobustGnl2, rsGnl2, Fmax] (intlinprog refers to max (r,s) robustness)
%

itertotal = length(testedMatrices);

%par
parfor ii=1:1:length(testedMatrices)
    
    TM = testedMatrices(ii,:);
    
    L = TM{1};
    %     graphnum = testedMatrices{ii,2};
    graphtype = TM{2};
    n = TM{3};
    pk = TM{4};
    algnum = TM{5};
    
    % LB algorithm (Algorithm 1 in Usevitch 2019 ACC)
    if algnum == 1
%         A = diag(diag(L)) - L;
        LBstruct = DetermineRobustness(struct('L',L,'smax',1));
        
        % Output ordering:
        % time r n p s
        
        LB = [LB; {LBstruct.time LBstruct.r n pk LBstruct.s graphtype}];
    end
    
    % DetRob
    if algnum == 2
        A = diag(diag(L)) - L;
        DetRobstruct = DetermineRobustness(struct('A',A));
        
        % Output ordering:
        % time r n p/k s
        
        DetRob = [DetRob; {DetRobstruct.time DetRobstruct.r n pk DetRobstruct.s graphtype}];
    end
    
    % rRobustGnl
    if algnum == 3
        rRGstruct = rRobustGnl(struct('L',L,'print',print));
        
        % Output ordering:
        % time r n p/k
        
        rRG = [rRG; {rRGstruct.time rRGstruct.maxr n pk graphtype}];
    end
    
    % rRobustGnl2
    if algnum == 4
        rRG2struct = rRobustGnl2(struct('L',L,'print',print));
        
        % Output ordering:
        % time r n p/k
        
        rRG2 = [rRG2; {rRG2struct.time rRG2struct.r n pk graphtype}];
    end
    
    % rsGnl2
    if algnum == 5
        rs2struct = rsGnl2(struct('L',L,'print',print));
        
        % Output ordering:
        % time r n p/k s
        
        rs2 = [rs2; {rs2struct.time rs2struct.r n pk rs2struct.s graphtype}];
    end
    
    % Fmax
    if algnum == 6
        Fmaxstruct = maxF(struct('L',L,'print',print));
        
        % Output ordering:
        % time r n p/k s
        % NOTE: F = r = s for this situation
        
        Fmax = [Fmax; {Fmaxstruct.time Fmaxstruct.Fhat n pk Fmaxstruct.Fhat graphtype}];
    end
    
    % rRoblower
    if algnum == 7
        rlowstruct = rRoblower(struct('L',L,'print',print));
        
        % Output ordering:
        % time r n p/k
        
        rlow = [rlow; {rlowstruct.time rlowstruct.r n pk graphtype}];
    end
    
    % rRobupper
    if algnum == 8
        rupstruct = rRobupper(struct('L',L,'print',print));
        
        % Output ordering:
        % time r n p/k
        
        rup = [rup; {rupstruct.time rupstruct.r n pk graphtype}];
    end
    
    if ii >= 10
        disp(['            Progress: row ' num2str(ii) ' completed, ' num2str(itertotal) ' total'])
    end
    
end

outstruct.testedMatrices = testedMatrices;
outstruct.LB = cell2mat(LB);
outstruct.DetRob = cell2mat(DetRob);
outstruct.rRG = cell2mat(rRG);
outstruct.rRG2 = cell2mat(rRG2);
outstruct.rs2 = cell2mat(rs2);
outstruct.Fmax = cell2mat(Fmax);
outstruct.rlow = cell2mat(rlow);
outstruct.rup = cell2mat(rup);

% outstruct.LB = LB;
% outstruct.DetRob = DetRob
% outstruct.rRG = rRG;
% outstruct.rRG2 = rRG2;
% outstruct.rs2 = rs2;
% outstruct.Fmax = Fmax;

end