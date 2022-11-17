function outstruct = rsGnl(args)

% General integer linear program (r,s)-robustness determination function.
% Determine the (r,s)-robustness of a digraph using intlinprog. Also measures
%   the computation time.
% args must be a struct with the following fields:
%   args.L : Graph Laplacian
%   args.print : (Optional) Controls whether output is displayed or not. Set
%       to 1 to print output of intlinprog and final results to the screen. Set
%       to 0 for no printed output. Default value is 0.
%   args.init : initial feasible point for the optimization problem

% NOTE!!!
% This function uses the ~(AvBvC) method to determine r.
% This function performs !! 2 !! optimization routines.
% The first determines (rho(D)+1,0), which determines rho(D).
% The second determines the max s for which the graph is (rho(D),s)-robust


L = args.L;

if size(L,1) ~= size(L,2)
    error('Matrix is not square')
end

n = size(L,1);
nr2 = ceil(n/2);

% Constraint matrices
AC = [0 0 -ones(1,n) zeros(1,3*n);
    0 0 ones(1,n) zeros(1,3*n);
    0 0 zeros(1,n) -ones(1,n) zeros(1,2*n);
    0 0 zeros(1,n) ones(1,n) zeros(1,2*n);
    zeros(n,1) zeros(n,1) eye(n) eye(n) zeros(n,n) zeros(n,n)];

BC = [-ones(n,1) zeros(n,1) L zeros(n,3*n);
    -ones(n,1) zeros(n,1) zeros(n) L zeros(n,2*n)];

CC = [0 0 -ones(1,n) zeros(1,n) ones(1,n) zeros(1,n);
    0 0 zeros(1,n) -ones(1,n) zeros(1,n) ones(1,n);
    0 -1 zeros(1,n) zeros(1,n) ones(1,n) ones(1,n)];

DC = [-ones(n,1) zeros(n,1) L zeros(n) -nr2*eye(n) zeros(n);
    -ones(n,1) zeros(n,1) zeros(n) L  zeros(n) -nr2*eye(n)];

EX = [0 0 zeros(1,2*n) -ones(1,2*n)]; % testing only

% A = [AC; BC; CC; DC];
A = [AC; CC; DC]; % testing only

bAC = [-1; (n-1); -1; (n-1); ones(n,1)];
bBC = [zeros(2*n,1)];
bCC = [-1;-1;0];
bDC = [-ones(2*n,1)];

% b = [bAC; bBC; bCC; bDC];
b = [bAC; bCC; bDC]; % Testing only
% b = [-1; (n-1); -1; (n-1); ones(n,1); zeros(2*n,1); -1; -1; 0; -ones(2*n,1)]; % orig;
% b = [-1; (n-1); -1; (n-1); ones(n,1); 2*(n-1)*ones(2*n,1); zeros(2*n,1); -1; -1; 0; -ones(n,1); -ones(n,1); -1]; % Testing


lb = [0; 0; zeros(4*n,1)];
ub = [Inf; Inf; ones(4*n,1)];

% c = [n+1; 1; zeros(4*n,1)];
c = [1; n+1; zeros(4*n,1)]; % testing only -- THIS WORKS

% Set Integer tolerance to the lowest possible value
% Change the Heuristics to basic or intermediate to possibly save time
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n),'Heuristics','advanced','MaxTime',1e20);

time = 0;
time2 = 0;

% Display or Suppress iterations, based on args.print
if ~isfield(args,'print') || args.print == 0
    options.Display = 'off';
end

if ~isfield(args,'init')
    tic
    [x,fval,exitflag,output] = intlinprog(c,3:4*n+2,A,b,[],[],lb,ub,options);
    time = toc;
else
    tic
    [x,fval,exitflag,output] = intlinprog(c,3:4*n+2,A,b,[],[],lb,ub,args.init,options);
    time = toc;
end

secondsolve = 0;
if x(2) == 0
    disp('Solving second problem for (r,s)')
    ub = [x(1)-1; Inf; ones(4*n,1)];
    
    tic
    [x2,fval2,exitflag2,output2] = intlinprog(c,3:4*n+2,A,b,[],[],lb,ub,options);
    time2 = toc;
    secondsolve = 1;
end

outstruct.b1 = x(3:n+2);
outstruct.b2 = x(n+3:2*n+2);
outstruct.y1 = x(2*n+3:3*n+2);
outstruct.y2 = x(3*n+3:end);
outstruct.exitflag = exitflag;
outstruct.output = output;
outstruct.time = time;

if secondsolve ==1
    outstruct.second.b1 = x2(3:n+2);
    outstruct.second.b2 = x2(n+3:2*n+2);
    outstruct.second.y1 = x2(2*n+3:3*n+2);
    outstruct.second.y2 = x2(3*n+3:end);
    outstruct.second.exitflag = exitflag2;
    outstruct.second.output = output2;
    outstruct.second.time = time2;
end

if secondsolve == 0
    
    outstruct.maxr = x(1);
    outstruct.maxs = x(2);
    
    if isfield(args,'print') && args.print == 1
        disp('Max r: ')
        x(1)
        disp('Max s: ')
        x(2)
        disp('Computation time:')
        time
        
        disp('One possible minimal S1 / S2 pair:')
        disp('S1:')
        outstruct.b1'
        disp('S2: ')
        outstruct.b2'
        
    end
else
    outstruct.maxr = x2(1);
    outstruct.maxs = x2(2);
    
    if isfield(args,'print') && args.print == 1
        disp('Max r: ')
        x2(1)
        disp('Max s: ')
        x2(2)
        disp('Computation time:')
        time + time2
        
        disp('One possible minimal S1 / S2 pair:')
        disp('S1:')
        outstruct.second.b1'
        disp('S2: ')
        outstruct.second.b2'
end




end