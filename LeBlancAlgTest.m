% Code to test LeBlanc algorithms

clear all
clc

args1 = struct;
args1.n = 20;
args1.k = 8;
args1.p = .8;
args1.type = 'kundir';

L = makegraph(args1);
A = -L + diag(diag(L));

args2 = struct('A',A,'smax',1);

tic
out = DetermineRobustness(args2);
time = toc

disp(strcat({'Time is: '},{' '},{num2str(time)}))


disp('R-Robustness:')
out.r

% Compare to ILP
n = args1.n;
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
optx = intlinprog(c',[2:n+1],A,-b);
toc

disp('Max r-Robustness:')
optx(1)

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
optx = intlinprog(c',[2:n+1],A,-b);
toc

disp('Max r-Robustness (int lin prog):')
optx(1)



