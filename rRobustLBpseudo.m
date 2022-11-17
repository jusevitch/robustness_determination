% r-Robustness Lower Bound
% This is an attempt to lower bound the r-Robustness of any arbitrary graph
% by using the inequality constraints from the equivalent LP.
%
% Letting x be the vector [t; sigma], we have
% Ax <= -b <=> -Ax >=b
% Since A has full column rank, we can use the Moore-Penrose pseudoinverse
% to obtain
% x >= -(A'*A)^(-1)*A'*b

clear all
clc

n = 45;
k = 20;

L = makegraph('kplatoon',n,k);

c = [1; zeros(n,1)];

% Ax + b <= 0 -- don't forget to make the b negative in intlinprog
% A = [-ones(n,1) -L;
%     -ones(n,1) L;
%     -1 zeros(1,n);
%     zeros(n,1) -eye(n);
%     zeros(n,1) eye(n);
%     0 -ones(1,n);
%     0 ones(1,n)];

% b = [zeros(n,1);
%     zeros(n,1);
%     0;
%     zeros(n,1);
%     -ones(n,1);
%     1
%     (1-n)];

% A,b with constraint to have 1st value 1
% A = [-ones(n,1) -L;
%     -ones(n,1) L;
%     -1 zeros(1,n);
%     zeros(n,1) -eye(n);
%     zeros(n,1) eye(n);
%     0 -ones(1,n);
%     0 ones(1,n);
%     0 c(1:end-1)'];
% 
% b = [zeros(n,1);
%     zeros(n,1);
%     0;
%     zeros(n,1);
%     -ones(n,1);
%     1
%     (1-n);
%     -1];


xbound = pinv(-A)*b



function outmatrix = makegraph(string,n,k)

if strcmp(string,'complete')
    %     Complete graph
    D = Dmatrix(n,[],'arbcomplete',0);
    outmatrix = D*D';
elseif strcmp(string,'kdir')
    % k-Circulant directed graph
    outmatrix = kCirculant(n,k,'dir');
elseif strcmp(string,'kundir')
    % k-Circulant undirected graph
    outmatrix = kCirculant(n,k,'undir');
elseif strcmp(string,'kplatoon')
    % k-nearest neighbor platoons;
    outmatrix = zeros(n);
    for ii=1:1:k
        outmatrix = outmatrix + diag(-ones(n - ii,1),ii) + diag(-ones(n-ii,1),-ii);
    end
    outmatrix = outmatrix + abs(diag(outmatrix*ones(n,1)));
else
    error('Sorry -- makegraph does not have that option')
end

end