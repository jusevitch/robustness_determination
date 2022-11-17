% rRobustWarmStartGurobi


clear all
clc

n = 100;
k = 49;
p = .8;

args = struct('n',n,'k',k,'p',p,'type','kundir');

L = makegraph(args);
dmin = min(diag(L));

display_final_sets = false; % Set to true if you want the final limiting sets displayed

%% Two-set case

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
ub = [dmin; ones(n,1)]; % Try both ub's to see which runs faster
% ub = [Inf; ones(2*n,1)];

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

tic
result = gurobi(model,params);
time1 = toc

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

disp('Incumbent robustness:')
fval

incsigma = incx(2:end);

if display_final_sets
    disp('Incumbent sigma value')
    incsigma'
end



%% General Case

clear A
clear b

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

% ***ADD THE FOLLOWING CONSTRAINTS:
% 0 \leq [I I]sigma \leq 1
% 0 <= t <= minimum in-degree of nodes in D

b = [(n-1)*ones(2*n,1); -(n-1)*ones(2*n,1); -1; (n-1); -1; (n-1); ones(n,1)];

lb = zeros(2*n+1,1);
% ub = [dmin + (n-1); ones(2*n,1)]; % Try both ub's to see which runs faster
ub = [Inf; ones(2*n,1)];

c = [1; zeros(2*n,1)];

% Set Integer tolerance to the lowest possible value
% options = optimoptions('intlinprog','IntegerTolerance',1e-6);

% [x,fval,exitflag,output] = intlinprog(c,1:2*n+1,A,b,[],[],lb,ub,options);

% Gurobi implementation

model.obj = c;
model.A = sparse(A);
model.vtype = repmat('I',2*n+1,1);
model.vtype(1) = 'C'; % for the variable t; this could also be an integer
model.sense = repmat('<',size(A,1),1); % Inequality constraints
model.rhs = full(b);
model.lb = lb;
model.ub = ub;

% Set the initial starting point
model.start = [incx(1) + (n-1); incsigma; (ones(n,1) - incsigma)]; % Remember to include the +(n-1) for the first entry

params = struct();
params.LogFile = 'gurobi_log.txt';
% params.LogToConsole = 0; % Suppresses screen output

tic
result = gurobi(model,params);
time2 = toc

% Resolve the model if status is infinite or unbounded
if strcmp(result.status,'INF_OR_UNBD')
    params.DualReductions = 0;
    warning('Infeasible or unbounded, resolve without dual reductions to determine...');
    result = gurobi(model,params);
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
    %         L = result.objbound;
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

disp('r-Robustness:')
fval - (n-1)
disp('Total optimization time:')
total_time = time1 + time2

if display_final_sets
    disp('One possible minimal S1 / S2 pair:')
    sigma = x(2:end);
    sigma1 = sigma(1:n)'
    sigma2 = sigma(n+1:end)'
end