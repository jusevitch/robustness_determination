% ILP for determining r-Robustness with nonempty S0 sets
% Theory and Algorithm developed by James Usevitch

clear all
clc

n = 50;
k = 20;

args = struct('n',n,'k',k,'type','kundir');

% L = makegraph('kundir',n,k);
L = makegraph(args);

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
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n));

[x,fval,exitflag,output] = intlinprog(c,1:2*n+1,A,b,[],[],lb,ub,options);

disp('Robustness:')
fval - (n-1)

sigma = x(2:end);
sigma1 = sigma(1:n);
sigma2 = sigma(n+1:end);

disp('One possible minimal S1 / S2 pair:')
disp('S1:')
sigma1'
disp('S2: ')
sigma2'


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

% A = [-ones(n,1) -L zeros(n,n);
%     -ones(n,1) zeros(n,n) -L;
%     -ones(n,1) L zeros(n,n);
%     -ones(n,1) zeros(n,n) L;
%     0 -ones(1,n) zeros(1,n);
%     0 ones(1,n) zeros(1,n);
%     0 zeros(1,n) -ones(1,n);
%     0 zeros(1,n) ones(1,n);
%     zeros(n,1) -eye(n) zeros(n);
%     zeros(n,1) zeros(n) -eye(n);
%     zeros(n,1) eye(n) zeros(n);
%     zeros(n,1) zeros(n) eye(n)];

% b = [(n-1)*ones(2*n,1); -(n-1)*ones(2*n,1); -1; (n-1); -1; (n-1); zeros(2*n,1); ones(2*n,1)];