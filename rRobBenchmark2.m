function outstruct = rRobBenchmark2(args)

% NOTE: Run this function at startup with a small nvec containing small values of n so
%   that Matlab can cache a parsed version. This will make it run faster
%   for subsequent function calls.
%   (see https://www.mathworks.com/matlabcentral/answers/415728-details-on-why-functions-are-faster-than-scripts)
%
% Argument fields:
% args.nvec : vector of n values to test. Column vector.
% args.pvec : vector of edge probability values to test. Row vector.
% args.poolsize : number of workers to include in the parpool
% args.graphs : Set entry to 1 to test the graph type, 0 otherwise. Graph types are (see makegraph.m function for details):
%               [Erdos, complete, kundir, kdir, randdir, kinrand, koutrand]
% args.algs : Set entry to 1 to test the algorithm group, 0 otherwise. Algorithm groups are [LB; intlinprog and intlinprogwm; Gurobi and Gurobiwm]
% args.testedMatrices : Matrix of pre-existing testedMatrices info. USE THIS WHEN CALLING THE GUROBI ALGORITHMS ONLY.
% args.testWarm : (GUROBI ONLY) Set to 0 if you do not want to test the
%                 Gurobi warm start algorithm
% args.test1robust : Ensure the graphs are at least 1 robust before testing
%                    them (checks Fiedler eigenvalue for undirected graphs,
%                    uses Tutte spanning tree test for digraphs)
% args.saveMatrices : Set to 1 to save the Laplacian matrices of the graphs tested in
%                     the variable outstruct.testedMatrices

finalLB = struct();
finalintlin = struct();
finalintlinwm = struct();
finalgurobi = struct();
finalgurobiwm = struct();

testedMatrices = [];

% Seed the random number generator to get different graphs
rng('shuffle');

if ~isfield(args,'graphs')
    graphs = [1;1;1;1;1;0;0]; % Test all graph types except for kinrand and koutrand
else
    graphs = args.graphs;
end

if ~isfield(args,'algs')
    algs = [1;1;1]; % Test all algorithms
else
    algs = args.algs;
end

if algs(1) || algs(2)
    
    graphslength = length(graphs);
    
    if graphs(1)
        Erdos = rRobBenchErdos(args);
        finalLB.Erdos = Erdos.LB;
        finalintlin.Erdos = Erdos.intlin;
        finalintlinwm.Erdos = Erdos.intlinwm;
        testedMatrices = [testedMatrices; Erdos.testedMatrices];
    end
    
    if graphs(2)
        complete = rRobBenchComplete(args);
        finalLB.complete = complete.LB;
        finalintlin.complete = complete.intlin;
        finalintlinwm.complete = complete.intlinwm;
        testedMatrices = [testedMatrices; complete.testedMatrices];
    end
    
    if graphs(3)
        kundir = rRobBenchkundir(args);
        finalLB.kundir = kundir.LB;
        finalintlin.kundir = kundir.intlin;
        finalintlinwm.kundir = kundir.intlinwm;
        testedMatrices = [testedMatrices; kundir.testedMatrices];
    end
    
    if graphs(4)
        kdir = rRobBenchkdir(args);
        finalLB.kdir = kdir.LB;
        finalintlin.kdir = kdir.intlin;
        finalintlinwm.kdir = kdir.intlinwm;
        testedMatrices = [testedMatrices; kdir.testedMatrices];
    end
    
    if graphslength > 4 && graphs(5) % Random directed graphs
        randdir = rRobBenchranddir(args);
        finalLB.randdir = randdir.LB;
        finalintlin.randdir = randdir.intlin;
        finalintlinwm.randdir = randdir.intlinwm;
        testedMatrices = [testedMatrices; randdir.testedMatrices];
    end
    
    if graphslength > 4 && graphs(6) % k-in random directed graphs
        kinrand = rRobBenchkinrand(args);
        finalLB.kinrand = kinrand.LB;
        finalintlin.kinrand = kinrand.intlin;
        finalintlinwm.kinrand = kinrand.intlinwm;
        testedMatrices = [testedMatrices; kinrand.testedMatrices];
    end
    
    if graphslength > 4 && graphs(7) % k-out random directed graphs
        koutrand = rRobBenchkoutrand(args);
        finalLB.koutrand = koutrand.LB;
        finalintlin.koutrand = koutrand.intlin;
        finalintlinwm.koutrand = koutrand.intlinwm;
        testedMatrices = [testedMatrices; koutrand.testedMatrices];
    end
    
    outstruct.finalLB = finalLB;
    outstruct.finalintlin = finalintlin;
    outstruct.finalintlinwm = finalintlinwm;
    
    %     testedMatrices = [Erdos.testedMatrices; complete.testedMatrices; kundir.testedMatrices; kdir.testedMatrices];
end

% if isempty(testedMatrices)
%     testedMatrices = args.testedMatrices;
% end

if (~algs(1) && ~algs(2)) && algs(3) % Only perform this code if the Gurobi algorithm is solely specified
    if ~isfield(args,'testedMatrices') || isempty(args.testedMatrices)
        % Create set of graphs to test with gurobiGraphs function
        testedMatrices = gurobiGraphs(args);
    else
        testedMatrices = args.testedMatrices;
    end
end

if algs(3)
    if isfield(args,'testWarm') && args.testWarm == 0
        argsg = struct('testWarm',args.testWarm);
    end
    %     if args.testWarm == 1 % Only done this way because of a stupid MATLAB bug that won't let me directly assign this.
    %         argsg.testWarm = 1;
    %     elseif args.testWarm == 0
    %         argsg.testWarm = 0;
    %     end
    argsg.testedMatrices = testedMatrices;
    gurobi = rRobBenchGurobi(argsg);
    outstruct.finalgurobi = gurobi.finalgurobi;
    outstruct.finalgurobiwm = gurobi.finalgurobiwm;
    
    if isfield(args,'saveMatrices') && args.saveMatrices == 1
        outstruct.testedMatrices = testedMatrices; % Passes out the cell matrix containing the tested matrices
    end
    
end


end