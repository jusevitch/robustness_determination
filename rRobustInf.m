function outstruct = rRobustInf(args)


% args must be a struct with the following fields:
%   args.L : Graph Laplacian
%   args.print : (Optional) Controls whether output is displayed or not. Set
%       to 1 to print output of intlinprog and final results to the screen. Set
%       to 0 for no printed output. Default value is 0.
%   args.init : initial feasible point for the optimization problem

L = args.L;

if size(L,1) ~= size(L,2)
    error('Matrix is not square')
end

n = size(L,1);
% nr2 = ceil(n/2);

A = [-ones(n,1) -L;
    -ones(n,1) L;
%     -1 zeros(1,n); % Unnecessary with lb constraints
%     zeros(n,1) -eye(n);
%     zeros(n,1) eye(n);
    0 -ones(1,n);
    0 ones(1,n)];

b = [zeros(n,1);
    zeros(n,1);
%     0;
%     zeros(n,1);
%     ones(n,1);
    -1
    (n-1)];

c = [1; zeros(n,1)];

lb = [0; zeros(n,1)];
ub = [Inf; ones(n,1)]; % Seems to run faster with the upper bound

% Set Integer tolerance to the lowest possible value
% Change the Heuristics to basic or intermediate to possibly save time
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n),'Heuristics','advanced','MaxTime',1e20);

time = 0;

% Display or Suppress iterations, based on args.print
if ~isfield(args,'print') || args.print == 0
    options.Display = 'off';
end


tic
[x,fval,exitflag,output] = intlinprog(c,2:n+1,A,b,[],[],lb,ub,options);
% [x,fval,exitflag,output] = linprog(c,[],[],Aeq,beq,lb,ub,options);
% [x,fval,exitflag,output] = linprog(c,A,b,[],[],lb,ub,options);
time = toc;


outstruct.r = x(1);
outstruct.b = x(2:n+1);
outstruct.exitflag = exitflag;
outstruct.output = output;
outstruct.time = time;

if isfield(args,'print') && args.print == 1
    disp(['Max r: ' num2str(x(1))])
    disp(['Computation time: ' num2str(time)])
    
    disp('One possible minimal S1 / S2 pair (S1 = ones(n,1) - S2):')
    round(outstruct.b')
end


end