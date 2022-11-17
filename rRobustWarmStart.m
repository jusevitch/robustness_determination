% r-Robustness ILP with Warm Start
% intlinprog version -- USE WITH MATLAB 2018a OR LATER

% This script solves the r-Robust problem in 2 steps:
%   1. Solve the ILP with S1 u S2 = V to obtain an incumbent value for the
%   general case
%   2. Solve the general case where S0 is possibly nonempty
% A good initial incumbent value can sometimes speed up the optimization
% process, i.e. fathoming nodes and eliminating the gap between (lower)
% best bound and the (upper bound) incumbent value.
%
% If cases where the r-robustness is determined when |S0| > 0 are rare
% (which may or may not be true), then this would end up being the right
% way to solve the problem anyways.
%
% If nothing else, we can more quickly obtain a tighter upper bound on the
% r-robustness of the graph.

clear all
clc

n = 30;
k = 8;
p = .4;

args = struct('n',n,'k',k,'p',p,'type','complete');

% L = makegraph('kundir',n,k);
L = makegraph(args);

c = [1; zeros(n,1)];

% Ax + b <= 0 -- don't forget to make the b negative in intlinprog
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

tic
incoptx = intlinprog(c',2:n+1,A,-b);
toc

disp('Incumbent r-robustness value:')
round(incoptx(1))

% Create incumbent point for general case
x0 = [incoptx; (ones(n,1) - incoptx(2:end))];
x0 = round(x0);
x0(1) = x0(1) + (n-1);

dmin = min(min(diag(L)), ceil(n/2));

% Constraint matrices
A = [-ones(n,1) -L zeros(n,n);
    -ones(n,1) zeros(n,n) -L;
    -ones(n,1) L zeros(n,n);
    -ones(n,1) zeros(n,n) L;
    zeros(n,1) -eye(n,n) -eye(n,n); % Guarantees 0 \leq [I I]sig \leq 1
    zeros(n,1) eye(n,n) eye(n,n); % See above
    0 -ones(1,n) zeros(1,n);
    0 ones(1,n) zeros(1,n);
    0 zeros(1,n) -ones(1,n);
    0 zeros(1,n) ones(1,n);
    zeros(n,1) eye(n) eye(n)];

% ADD THE CONSTRAINT THAT 0 <= t <= minimum in-degree of nodes in D
% This is because the r-robustness of the graph is upper bounded by the
% minimum in-degree of the nodes

b = [(n-1)*ones(2*n,1); -(n-1)*ones(2*n,1); zeros(n,1); ones(n,1); -1; (n-1); -1; (n-1); ones(n,1)];

lb = zeros(2*n+1,1);
ub = [Inf; ones(2*n,1)]; % TESTING PURPOSES ONLY
% ub = [dmin + (n-1); ones(2*n,1)]; % With upper bound of t <= minimum in-degree of nodes in D
% The (n-1) needs to be added since the objective minimizes t + (n-1)
% The program seems to run faster when the upper bound on t is excluded.

c = [1; zeros(2*n,1)];

% Set Integer tolerance to the lowest possible value
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n),'Heuristics','advanced','MaxTime',1e20);

tic
[x,fval,exitflag,output] = intlinprog(c,2:2*n+1,A,b,[],[],lb,ub,x0,options);
toc

disp('r-Robustness:')
fval - (n-1)

sigma = x(2:end);
sigma1 = sigma(1:n);
sigma2 = sigma(n+1:end);

disp('One possible minimal S1 / S2 pair:')
disp('S1:')
sigma1'
disp('S2: ')
sigma2'