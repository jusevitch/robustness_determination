function outstruct = maxF(args)

% Determines the maximum value of F for which a digraph is (F+1,F+1)-robust. Also measures
%   the computation time.
% args must be a struct with the following fields:
%   args.L : Graph Laplacian
%   args.print : (Optional) Controls whether output is displayed or not. Set
%       to 1 to print output of intlinprog and final results to the screen. Set
%       to 0 for no printed output. Default value is 0.
%   args.init : initial feasible point for the optimization problem
%   args.MaxTime : Maximum time limit for algorithm to run. Do not define
%                   this field for default value of 1e20


L = args.L;

if size(L,1) ~= size(L,2)
    error('Matrix is not square')
end

time = 0;

tic

n = size(L,1);
nr2 = ceil(n/2);

AC = [0 -ones(1,n) zeros(1,3*n);
    0 ones(1,n) zeros(1,3*n);
    0 zeros(1,n) -ones(1,n) zeros(1,2*n);
    0 zeros(1,n) ones(1,n) zeros(1,2*n);
    zeros(n,1) eye(n) eye(n) zeros(n,2*n)];

BC = [0 -ones(1,n) zeros(1,n) ones(1,n) zeros(1,n);
    0 zeros(1,n) -ones(1,n) zeros(1,n) ones(1,n);
    -1 zeros(1,2*n) ones(1,2*n)];

CC = [-ones(n,1) L zeros(n) -nr2*eye(n) zeros(n);
    -ones(n,1) zeros(n) L zeros(n) -nr2*eye(n)];

DC = [zeros(n,1) zeros(n,2*n) eye(n) eye(n)]; % Constrains y1 + y2 <= 1

% A = [AC; BC; CC];
A = [AC; BC; CC; DC]; % Testing

bAC = [-1; n-1; -1; n-1; ones(n,1)];

bBC = [-1; -1; 0];

bCC = -zeros(2*n,1);

bDC = ones(n,1);

% b = [bAC; bBC; bCC];
b = [bAC; bBC; bCC; bDC]; % Testing

lb = [-1; zeros(4*n,1)];
ub = [ceil(n/2)-1; ones(4*n,1)];

c = [1; zeros(4*n,1)];

% Set Integer tolerance to the lowest possible value
% Change the Heuristics to basic or intermediate to possibly save time
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n),'Heuristics','advanced');
if isfield(args,'MaxTime')
    options.MaxTime = args.MaxTime;
else
    options.MaxTime = 1e20;
end



% Display or Suppress iterations, based on args.print
if ~isfield(args,'print') || args.print == 0
    options.Display = 'off';
end


[x,fval,exitflag,output] = intlinprog(c,2:4*n+1,A,b,[],[],lb,ub,options);
% [x,fval,exitflag,output] = linprog(c,[],[],Aeq,beq,lb,ub,options);
% [x,fval,exitflag,output] = linprog(c,A,b,[],[],lb,ub,options);
time = toc;

% Note: the optimization routine finds the min value of Fbar such that the
% graph is NOT (Fbar,Fbar)-robust. Therefore the max value of p such
% that the graph is (p,p)-robust is p = Fbar-1.
% if isempty(x)
%     disp('here')
% end
% *** THE OPTIMIZATION PROBLEM SEEMS TO BE INFEASIBLE FOR COMPLETE GRAPHS

if ~isempty(x)
    outstruct.Fhat = x(1);
    if x(1) > 0
        outstruct.F = x(1)-1;
    else
        outstruct.F = 0;
    end
    outstruct.b1 = x(2:n+1);
    outstruct.b2 = x(n+2:2*n+1);
    outstruct.y1 = x(2*n+2:3*n+1);
    outstruct.y2 = x(3*n+2,end);
    outstruct.exitflag = exitflag;
    outstruct.output = output;
    outstruct.time = time;
    
    if isfield(args,'print') && args.print == 1
        disp(['Max value of Fhat: ' num2str(x(1))])
        disp(['Max value of F: ' num2str(x(1)-1)])
        disp(['Computation time: ' num2str(time)])
    end
end

if exitflag ~= 1
    disp(['Solution did not converge--exitflag = ' num2str(exitflag)]);
end

end