function outstruct = rRobustGnl(args)

% General integer linear program r-Robustness determination function.
% Determine the r-Robustness of a digraph using intlinprog. Also measures
%   the computation time.
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
% Change the Heuristics to basic or intermediate to possibly save time
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n),'Heuristics','advanced','MaxTime',1e20);

% Display or Suppress iterations, based on args.print
if ~isfield(args,'print') || args.print == 0
    options.Display = 'off';
end

tic
[x,fval,exitflag,output] = intlinprog(c,2:2*n+1,A,b,[],[],lb,ub,options);
time = toc;

sigma = x(2:end);

outstruct.maxr = fval - (n-1);
outstruct.sigma1 = x(2:n+1);
outstruct.sigma2 = x(n+2:end);
outstruct.exitflag = exitflag;
outstruct.output = output;
outstruct.time = time;

if isfield(args,'print') && args.print == 1
    disp('Robustness:')
    fval - (n-1)
    disp('Computation time:')
    time
    
    disp('One possible minimal S1 / S2 pair:')
    disp('S1:')
    outstruct.sigma1'
    disp('S2: ')
    outstruct.sigma2'
    
end

end