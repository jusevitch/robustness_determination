% Gurobi Implementation of the maximum r-robustness determination problem
% Theory by James Usevitch; code by James Usevitch and the Gurobi Tutorial
% Website

clear all
clc

n = 21;
k = 20;
p = .8;

args = struct('n',n,'k',k,'p',p,'type','complete');

L = makegraph(args);

% Random matrix methodology
% Entries of Adjacency matrix chosen randomly
% A better method should be implemented that makes typical random graphs such as
% Barabasi-Albert, Erdos-Renyi, etc.
% Adj = randi([0,1],n,n);
% D = diag(Adj*ones(n,1));
% L = D - Adj;

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

params = struct();
params.LogFile = 'gurobi_log.txt';
% params.LogToConsole = 0; % Suppresses screen output

result = gurobi(model,params);

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

disp('Robustness:')
fval - (n-1)

sigma = x(2:end);
sigma1 = sigma(1:n);
sigma2 = sigma(n+1:end);


% function outmatrix = makegraph(string,n,k)
% 
% if strcmp(string,'complete')
%     %     Complete graph
%     D = Dmatrix(n,[],'arbcomplete',0);
%     outmatrix = D*D';
% elseif strcmp(string,'kdir')
%     % k-Circulant directed graph
%     outmatrix = kCirculant(n,k,'dir');
% elseif strcmp(string,'kundir')
%     % k-Circulant undirected graph
%     outmatrix = kCirculant(n,k,'undir');
% elseif strcmp(string,'kplatoon')
%     % k-nearest neighbor platoons;
%     outmatrix = zeros(n);
%     for ii=1:1:k
%         outmatrix = outmatrix + diag(-ones(n - ii,1),ii) + diag(-ones(n-ii,1),-ii);
%     end
%     outmatrix = outmatrix + abs(diag(outmatrix*ones(n,1)));
% else
%     error('Sorry -- makegraph does not have that option')
% end
% 
% end