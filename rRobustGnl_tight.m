function outstruct = rRobustGnl_tight(args)

% Integer linear program r-Robustness determination function with tighter constraints.
% Determine the r-Robustness of a digraph using intlinprog. Also measures
%   the computation time.
% args must be a struct with the following fields:
%   args.L : Graph Laplacian
%   args.print : (Optional) Controls whether output is displayed or not. Set
%       to 1 to print output of intlinprog and final results to the screen. Set
%       to 0 for no printed output. Default value is 0.

% NOTE:
% THIS METHOD SEEMS TO BE MUCH SLOWER THAN THE GNL METHOD ACCORDING TO
% INITIAL TESTS. Try the rRobustGnl formulation as well to compare speeds.

L = args.L;

if size(L,1) ~= size(L,2)
    error('Matrix is not square')
end

n = size(L,1);

% *** DANGER *** 
% FOR DIRECTED GRAPHS, ONE NODE CAN HAVE ZERO IN-DEGREE, BUT THE GRAPH IS
% STILL 1-ROBUST
% CONSIDER A ROOTED OUTBRANCHING
dmin = min(max(min(diag(L)),1), ceil(n/2));

% Vector of powers of 2 - SCALES BADLY
% twovec = 0:1:n-1;
% twovec = twovec';
% twovec = 2.^twovec;

% Matrix of powers of two
% Tmat = diag(twovec);

% Enumeration vector
numvec = 1:n;
numvec = numvec';

% numvec repeated
NM = kron(numvec',ones(n,1));

% Variable vector: [t s sigma_0' sigma_1']'
% Size is 2n + 2

% Constraint matrices
A = [-ones(n,1) zeros(n,n) zeros(n,n) -L ; % constraint A1
    -ones(n,1) zeros(n,n) zeros(n,n) L; % A2
    -ones(n,1) zeros(n,n) L L; % A3
    -ones(n,1) zeros(n,n) -L -L; % %A4 end part A
    zeros(n,1) zeros(n,n) -eye(n,n) zeros(n,n); % B1 begin part B
    zeros(n,1) zeros(n,n) eye(n,n) zeros(n,n); % B2
    zeros(n,1) zeros(n,n) zeros(n,n) -eye(n,n); % B3
    zeros(n,1) zeros(n,n) zeros(n,n) eye(n,n); % B4
    0 zeros(1,n) ones(1,n) zeros(1,n); % constraint B6 (B5 was unnecessary)
    0 zeros(1,n) zeros(1,n) -ones(1,n); % B7
    0 zeros(1,n) zeros(1,n) ones(1,n); % B8
    0 zeros(1,n) ones(1,n) ones(1,n); % B9
    0 zeros(1,n) -ones(1,n) -ones(1,n); % B10
    zeros(n,1) zeros(n,n) -eye(n,n) -eye(n,n); % B11
    zeros(n,1) zeros(n,n) eye(n,n) eye(n,n); % B12 end part B
    zeros(n,1) -NM -diag(numvec) zeros(n,n); % C1 begin part C
    zeros(n,1) eye(n,n) zeros(n,n) -eye(n,n)]; % C2
    

% ADD THE CONSTRAINT THAT 0 <= t <= minimum in-degree of nodes in D?
% This is because the r-robustness of the graph is upper bounded by the
% minimum in-degree of the nodes

b = [(n-1)*ones(n,1); % constraint A1
    -(n-1)*ones(n,1); % A2
    (n-1)*ones(n,1); % A3
    -(n-1)*ones(n,1); % A4
    zeros(n,1); % B1
    ones(n,1); % B2
    zeros(n,1); % B3
    ones(n,1); % B4
    n-2; % constraint B6
    -1; % B7
    n-1; % B8
    n-1; % B9
    -1; % B10
    zeros(n,1); % B11
    ones(n,1); % B12 end part B
    -numvec; % C1 begin part C
    zeros(n,1)]; % C2

% Equality Constraint Matrix

Aeq = [0 ones(1,n) zeros(1,n) zeros(1,n)];

beq = 1;

lb = [0; zeros(n,1); zeros(2*n,1)];
ub = [Inf; ones(n,1); ones(2*n,1)];
% ub = [dmin + (n-1); ones(2*n,1)]; % With upper bound of t <= minimum in-degree of nodes in D
% The (n-1) needs to be added since the objective minimizes t + (n-1)
% The program seems to run faster when the upper bound on t is excluded.

c = [1; zeros(n,1); zeros(2*n,1)];

% Set Integer tolerance to the lowest possible value
% Change the Heuristics to basic or intermediate to possibly save time
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n),'Heuristics','advanced','MaxTime',1e20);

% Display or Suppress iterations, based on args.print
if ~isfield(args,'print') || args.print == 0
    options.Display = 'off';
end

tic
[x,fval,exitflag,output] = intlinprog(c,2:3*n+1,A,b,Aeq,beq,lb,ub,options);
time = toc;

sigma = x(2:end);

outstruct.maxr = fval - (n-1);
% Sigma values are different; see variable organization above
% sigma2 = ones(n,1) - sigma_0 - sigma_1
outstruct.sigma1 = x(2*n+2:end);
outstruct.sigma2 = ones(n,1) - x(n+2:2*n+1) - outstruct.sigma1;
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

