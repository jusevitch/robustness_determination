function outstruct = rRobustWm(args)

% Warm start integer linear program r-Robustness determination function.
% Determine the r-Robustness of a digraph using intlinprog and the two-set 
%   method as a warm start. Also measures the computation time.
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

% dmin = min(min(diag(L)), ceil(n/2));

%% Two-set Case

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

% Set Integer tolerance to the lowest possible value
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n));

% Display or Suppress iterations, based on args.print
if ~isfield(args,'print') || args.print == 0
    options.Display = 'off';
end

tic
incoptx = intlinprog(c',[2:n+1],A,-b,[],[],[],[],options);
time1 = toc;

if isfield(args,'print') && args.print == 1
disp('Incumbent r-robustness value:')
end

round(incoptx(1));

% Create incumbent point for general case
x0 = [incoptx; (ones(n,1) - incoptx(2:end))];
x0 = round(x0);
x0(1) = x0(1) + (n-1);

%% General Case

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

b = [(n-1)*ones(2*n,1); -(n-1)*ones(2*n,1); zeros(n,1); ones(n,1); -1; (n-1); -1; (n-1); ones(n,1)];

lb = zeros(2*n+1,1);
ub = [Inf; ones(2*n,1)]; % TESTING PURPOSES ONLY
% ub = [dmin + (n-1); ones(2*n,1)]; % With upper bound of t <= minimum in-degree of nodes in D
% The (n-1) needs to be added since the objective minimizes t + (n-1)
% The program seems to run faster when the upper bound on t is excluded.

c = [1; zeros(2*n,1)];

% WARNING: THIS WILL NOT WORK ON ANYTHING EARLIER THAN MATLAB 2018a
tic
[x,fval,exitflag,output] = intlinprog(c,2:2*n+1,A,b,[],[],lb,ub,x0,options); 
time2 = toc;

sigma = x(2:end);

outstruct.maxr = fval - (n-1);
outstruct.sigma1 = x(2:n+1);
outstruct.sigma2 = x(n+2:end);
outstruct.exitflag = exitflag;
outstruct.output = output;
outstruct.time = time1 + time2;
outstruct.time1 = time1;
outstruct.time2 = time2;

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