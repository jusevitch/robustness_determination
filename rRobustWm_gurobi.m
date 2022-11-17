function outstruct = rRobustWm_gurobi(args)

% Warm start integer linear program r-Robustness determination
% function using Gurobi.
% Determine the r-Robustness of a digraph using the Gurobi solver and the 
%   two-set warm start method. Also measures the computation time.
% THE GUROBI SOLVER MUST BE INSTALLED BEFORE RUNNING THIS FUNCTION.
%
% args must be a struct with the following fields:
%   args.L : Graph Laplacian
%   args.print : (Optional) Controls whether output is displayed or not. Set
%       to 1 to print output of intlinprog and final results to the screen. Set
%       to 0 for no printed output. Default value is 0.

L = args.L;

if size(L,1) ~= size(L,2)
    error('Matrix is not square')
end

n = size(L,1);

%% Two-set Case

% Constraint matrices
% Ax + b <= 0 -- you may need to make b negative if Gurobi requires Ax <= -b

A = [-ones(n,1) -L;
    -ones(n,1) L;
    -1 zeros(1,n);
    zeros(n,1) -eye(n);
    zeros(n,1) eye(n);
    0 -ones(1,n);
    0 ones(1,n)];

b = [zeros(n,1);
    zeros(n,1);
    0;
    zeros(n,1);
    -ones(n,1);
    1
    (1-n)];

lb = zeros(n+1,1);
% ub = [dmin; ones(n,1)]; % Try both ub's to see which runs faster
ub = [Inf; ones(n,1)];

c = [1; zeros(n,1)];

model.obj = c;
model.A = sparse(A);
model.vtype = repmat('I',n+1,1);
model.vtype(1) = 'C'; % for the variable t; this could also be an integer
model.sense = repmat('<',size(A,1),1); % Inequality constraints
model.rhs = full(-b);
model.lb = lb;
model.ub = ub;

params = struct();

% Suppress console output by default
if ~isfield(args,'print') || args.print == 0
    params.OutputFlag = 0;
end

tic
result = gurobi(model,params);
time1 = toc;

% Resolve the model if status is infinite or unbounded
if strcmp(result.status,'INF_OR_UNBD')
    params.DualReductions = 0;
    warning('Infeasible or unbounded, resolve without dual reductions to determine...');
    result = gurobi(model,params);
end

% Get the results
if isfield(result,'x')
    incx = result.x;
    %     if nargout > 3
    %         slack = model.A*x-model.rhs;
    %         violA = slack(1:size(A,1));
    %         violAeq = norm(slack((size(A,1)+1):end),inf);
    %         viollb = model.lb(:)-x;
    %         violub = 0;
    %         if isfield(model,'ub')
    %             violub = x-model.ub(:);
    %         end
    %         output.constrviolation = max([0; violA; violAeq; viollb; violub]);
    %     end
end

fval = [];

if isfield(result,'objval')
    fval = result.objval;
    %     if nargout > 3 && numel(intcon) > 0
    %         U = fval;
    %         L = result.objbound;
    %         output.relativegap = 100*(U-L)/(abs(U)+1);
    %         output.absolutegap = U-L;
    %     end
end

% if strcmp(result.status, 'OPTIMAL')
%     exitflag = 1;
% elseif strcmp(result.status, 'INFEASIBLE') ...
%         || strcmp(result.status, 'CUTOFF')
%     exitflag = -2;
% elseif strcmp(result.status, 'UNBOUNDED')
%     exitflag = -3;
% elseif isfield(result, 'x')
%     exitflag = 2;
% else
%     exitflag = 0;
% end

incsigma = incx(2:end);

if isfield(args,'print') && args.print == 1
disp('Incumbent robustness:')
fval

disp('Incumbent sigma value')
    incsigma'
end



%% General Case

% Constraint matrices
A = [-ones(n,1) -L zeros(n,n);
    -ones(n,1) zeros(n,n) -L;
    -ones(n,1) L zeros(n,n);
    -ones(n,1) zeros(n,n) L;
    0 -ones(1,n) zeros(1,n);
    0 ones(1,n) zeros(1,n);
    0 zeros(1,n) -ones(1,n);
    0 zeros(1,n) ones(1,n);
    zeros(n,1) eye(n) eye(n)];

b = [(n-1)*ones(2*n,1); -(n-1)*ones(2*n,1); -1; (n-1); -1; (n-1); ones(n,1)];

% Objective Function
c = [1; zeros(2*n,1)];

% Bounds
lb = zeros(2*n+1,1);
ub = [Inf; ones(2*n,1)]; % First ub can be changed to minimum in-degree in the graph (this is upper bound on robustness)

% Set Integer tolerance to the lowest possible value
options = optimoptions('intlinprog','IntegerTolerance',1e-6);

% Gurobi implementation

model.obj = c;
model.A = sparse(A);% Constraint matrices
model.vtype = repmat('I',2*n+1,1);
model.vtype(1) = 'C'; % for the variable t; this could also be an integer
model.sense = repmat('<',size(A,1),1); % Inequality constraints
model.rhs = full(b);
model.lb = lb;
model.ub = ub;

% Set the initial starting point
model.start = [incx(1) + (n-1); incsigma; (ones(n,1) - incsigma)]; % Remember to include the +(n-1) for the first entry

params = struct();

% Suppress console output by default
if ~isfield(args,'print') || args.print == 0
    params.OutputFlag = 0;
end

tic
result = gurobi(model,params);
time2 = toc;

% Resolve the model if status is infinite or unbounded
if strcmp(result.status,'INF_OR_UNBD')
    params.DualReductions = 0;
    warning('Infeasible or unbounded, resolve without dual reductions to determine...');
end

% Get the results
if isfield(result,'x')
    x = result.x;
    %     if nargout > 3
    %         slack = model.A*x-model.rhs;
    %         violA = slack(1:size(A,1));
    %         violAeq = norm(slack((size(A,1)+1):end),inf);
    %         viollb = model.lb(:)-x;
    %         violub = 0;
    %         if isfield(model,'ub')
    %             violub = x-model.ub(:);
    %         end
    %         output.constrviolation = max([0; violA; violAeq; viollb; violub]);
    %     end
end

fval = [];

if isfield(result,'objval')
    fval = result.objval;
    %     if nargout > 3 && numel(intcon) > 0
    %         U = fval;
    %         L = result.objbound;incsigma = incx(2:end);
    %         output.relativegap = 100*(U-L)/(abs(U)+1);
    %         output.absolutegap = U-L;
    %     end
end

if strcmp(result.status, 'OPTIMAL')
    exitflag = 1;
elseif strcmp(result.status, 'INFEASIBLE') ...
        || strcmp(result.status, 'CUTOFF')
    exitflag = -2;
elseif strcmp(result.status, 'UNBOUNDED')
    exitflag = -3;
elseif isfield(result, 'x')
    exitflag = 2;
else
    exitflag = 0;
end

outstruct.maxr = fval - (n-1);
outstruct.sigma1 = x(2:n+1);
outstruct.sigma2 = x(n+2:end);
outstruct.exitflag = exitflag;
% outstruct.output = output; % Can allow this later
outstruct.time = time1 + time2;
outstruct.time1 = time1;
outstruct.time2 = time2;
outstruct.incumbent = incsigma; % Solution of the two-set case

if isfield(args,'print') && args.print == 1
    disp('Robustness:')
    fval - (n-1)
    disp('Computation time:')
    time1 + time2
    
    disp('One possible minimal S1 / S2 pair:')
    disp('S1:')
    outstruct.sigma1'
    disp('S2: ')
    outstruct.sigma2'
end

end