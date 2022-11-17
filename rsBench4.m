function outstruct = rsBench4(args)

% This version of rsBench uses parfeval instead of parfor. The code
% distributes tasks between workers much more efficiently.

% To plot the results, use rsBenchplot2.m or later version.

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
% args.MaxTime : Maximum amount of time for algorithm to run. This only
%                   applies to the MILP formulations and not to
%                   DetRob or LB algorithms.

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
algnames = {'LB','DetRob','rRG','rRG2','rs2','Fmax','rlow','rup'};

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
if isfield(args,'pvec')
    if ~isempty(args.pvec)
        pvec = args.pvec;
        pvec_exists = true;
    else
        warning('Empty pvec input.')
        pvec_exists = false;
    end
else
    if graphs(1) || graphs(2)
        error('Analysis of Erdos or randdir graphs requested, but no pvec field specified.')
    else
        pvec_exists = false;
    end
end
if isfield(args,'kvec')
    if ~isempty(args.kvec)
        kvec = args.kvec;
        kvec_exists = true;
    else
        warning('Empty kvec input')
        kvec_exists = false;
    end
else
    if graphs(3) || graphs(4)
        error('Analysis of kin or kout graphs requested, but no kvec field specified.')
    else
        kvec_exists = false;
    end
end
avec = find(algs); % Number representation of algorithms

% Long vector of all n and p combinations; e.g. [n1 p1; n1 p2; ... n1 pend; n2 p1; ...]
% Makes it possible for the parfor loop to split each scenario between the
% workers more efficiently so no one core gets bogged down by calculating
% all graphs with a large n value or a large p value
if pvec_exists
    nones = [ones(length(pvec),1) zeros(length(pvec),1)];
    pones = [zeros(length(nvec),1) ones(length(nvec),1)];
    totvecp = kron(nvec,nones) + kron(pones,pvec); % Vector of all [n p] combos
end


% Long vector of all n and k combinations; e.g. [n1 k1; n1 k2; ... n1 kend; n2 k1; ...]
if kvec_exists
    nones = [ones(length(kvec),1) zeros(length(kvec),1)];
    kones = [zeros(length(nvec),1) ones(length(nvec),1)];
    totveck = kron(nvec,nones) + kron(kones,kvec); % Vector of all [n p] combos
end

% Create master list of total matrices

lastpidx = 2; % index of the last graph requiring p parameter in graphs vector

for ii=1:1:length(graphs)
    graphtype = graphnames(ii);
    if ii <= lastpidx
        if pvec_exists
            totvec = totvecp;
            totvec_exists = true;
        else
            totvec_exists = false;
        end
    else
        if kvec_exists
            totvec = totveck;
            totvec_exists = true;
        else
            totvec_exists = false;
        end
    end
    
    if totvec_exists
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
end

% Repeat ALL matrices once for each algorithm speficied by args.algs
testedMatrices = [repelem(testedMatrices,length(avec),1) num2cell(repmat(avec,length(testedMatrices),1))];

outstruct.testedMatrices = testedMatrices; % Perhaps make this in separate function?

%
% Perform Algorithm testing
% Algorithm groups are [LB, DetRob, rRobustGnl, rRobustGnl2, rsGnl2, Fmax] (intlinprog refers to max (r,s) robustness)
%

q = parallel.pool.DataQueue();
afterEach(q, @disp);

for jj=1:1:length(algnames)
    outstruct.(algnames{jj}) = zeros(sum(cell2mat(testedMatrices(:,5)) == jj),9); % Only 7 columns since we'll exclude the algnumber column
end

itertotal = length(testedMatrices);

h = waitbar(0,'Constructing F object...');

for ii=1:itertotal
    if isfield(args,'MaxTime')
        F(ii) = parfeval(pool,@evalInstance,1,testedMatrices(ii,:),q,args.MaxTime); % only q at end
    else
        F(ii) = parfeval(pool,@evalInstance,1,testedMatrices(ii,:),q); % only q at end
    end
    waitbar(ii/itertotal,h);
end

waitbar(0,h,'Waiting for problem instances to complete...');

for kk=1:itertotal
    [completeidx, nextvec] = fetchNext(F);
    % Add nextvec to the appropriate outstruct field
    outstruct.(algnames{nextvec(end)})(find(outstruct.(algnames{nextvec(end)}) == zeros(1,9),1,'first'),:) = nextvec(1:end-1);
    waitbar(kk/itertotal,h,'Processing. This could take a while.')
    disp([num2str(kk) ' / ' num2str(itertotal)]);
end

close all hidden

% outstruct.testedMatrices = testedMatrices;
% outstruct.LB = cell2mat(LB);
% outstruct.DetRob = cell2mat(DetRob);
% outstruct.rRG = cell2mat(rRG);
% outstruct.rRG2 = cell2mat(rRG2);
% outstruct.rs2 = cell2mat(rs2);
% outstruct.Fmax = cell2mat(Fmax);
% outstruct.rlow = cell2mat(rlow);
% outstruct.rup = cell2mat(rup);

end



function outvector = evalInstance(incellmatrix,q,varargin)

% Subfunction to evaluate each problem instance
% The variable incellmatrix should be one row of testedMatrices
% varargin should be the MaxTime parameter for the MILP formulations.
% THE MaxTime PARAMETER DOES NOT APPLY TO DetermineRobustness

L = incellmatrix{1};
%     graphnum = testedMatrices{ii,2};
graphtype = incellmatrix{2};
n = incellmatrix{3};
pk = incellmatrix{4};
algnum = incellmatrix{5};

is_MaxTime_set = false;

if nargin == 3
    is_MaxTime_set = true;
    MaxTime = varargin{1};
end

% send(q,sprintf('MaxTime for worker: %d',MaxTime));

% task = getCurrentTask;
% send(q, sprintf('Worker %d',id));

% Data format:
% [time, r, s, n, pk, graphtype, algtype]
% algtype:
%   1: LB
%   2: DetRob
%   3: rRG
%   4: rRG2
%   5: rs2
%   6: Fmax
%   7: rlow
%   8: rup
% graphtype:
%   1: Erdos
%   2: randdir
%   3: kin
%   4: kout

% Data Format:
% [time, r, s, n, pk, graphtype, success_before_MaxTime_bool, r_ub, s_ub algnum]

% successbool is 1 if exact robustness found within time limits; 0
%   otherwise
% If the algorithm did not find optimal value before MaxTime, then r and s
%   are the approx lower bounds on r and s via branch-and-bound, and r_ub, s_ub are approx upper
%   bounds. If r and s were found before MaxTime, then r_ub, s_ub are equal
%   to r and s.
% For DetRob and LB, r_ub = r and s_ub = s.

if is_MaxTime_set
    switch algnum
        case 1 % LB algorithm (Algorithm 1 in Usevitch 2019 ACC)
            %         tempstruct = DetermineRobustness(struct('L',L,'smax',1,'q',q));
            tempstruct = DetermineRobustness(struct('L',L,'smax',1));
            outvector = [tempstruct.time tempstruct.r tempstruct.s n pk graphtype 1 tempstruct.r tempstruct.s algnum];
        case 2 % DetRob
            %         tempstruct = DetermineRobustness(struct('L',L,'q',q));
            tempstruct = DetermineRobustness(struct('L',L));
            outvector = [tempstruct.time tempstruct.r tempstruct.s n pk graphtype 1 tempstruct.r tempstruct.s algnum];
        case 3 % rRobustGnl
            % !!!! PROBABLY DOES NOT WORK!!!!!
            tempstruct = rRobustGnl(struct('L',L));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
        case 4 % rRobustGnl2
            tempstruct = rRobustGnl2(struct('L',L,'MaxTime',MaxTime));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
        case 5 % rsGnl2
            tempstruct = rsGnl2(struct('L',L,'MaxTime',MaxTime));
            outvector = [tempstruct.totaltime tempstruct.r tempstruct.s n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub tempstruct.s_ub algnum];
        case 6 % Fmax
            % !!!! DOES NOT WORK!!!!!
            tempstruct = maxF(struct('L',L,'MaxTime',MaxTime));
            outvector = [tempstruct.time tempstruct.Fhat tempstruct.Fhat n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.Fhat tempstruct.Fhat algnum];
        case 7 % rRoblower
            tempstruct = rRoblower(struct('L',L,'MaxTime',MaxTime));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
        case 8 % rRobupper
            tempstruct = rRobupper(struct('L',L,'MaxTime',MaxTime));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
    end
    
else
    switch algnum
        case 1 % LB algorithm (Algorithm 1 in Usevitch 2019 ACC)
            %         tempstruct = DetermineRobustness(struct('L',L,'smax',1,'q',q));
            tempstruct = DetermineRobustness(struct('L',L,'smax',1));
            outvector = [tempstruct.time tempstruct.r tempstruct.s n pk graphtype 1 tempstruct.r tempstruct.s algnum];
        case 2 % DetRob
            %         tempstruct = DetermineRobustness(struct('L',L,'q',q));
            tempstruct = DetermineRobustness(struct('L',L));
            outvector = [tempstruct.time tempstruct.r tempstruct.s n pk graphtype 1 tempstruct.r tempstruct.s algnum];
        case 3 % rRobustGnl
            tempstruct = rRobustGnl(struct('L',L));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
        case 4 % rRobustGnl2
            tempstruct = rRobustGnl2(struct('L',L));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
        case 5 % rsGnl2
            tempstruct = rsGnl2(struct('L',L));
            outvector = [tempstruct.totaltime tempstruct.r tempstruct.s n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub tempstruct.s_ub algnum];
        case 6 % Fmax
            tempstruct = maxF(struct('L',L));
            outvector = [tempstruct.time tempstruct.Fhat tempstruct.Fhat n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.Fhat tempstruct.Fhat algnum];
        case 7 % rRoblower
            tempstruct = rRoblower(struct('L',L));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
        case 8 % rRobupper
            tempstruct = rRobupper(struct('L',L));
            outvector = [tempstruct.time tempstruct.r 1 n pk graphtype tempstruct.success_before_MaxTime_bool tempstruct.r_ub Inf algnum];
    end
end


end