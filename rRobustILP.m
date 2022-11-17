% r-Robust MATLAB Linear Programming

% THIS FILE TREATS THE CASE WHEN S1 u S2 = V
% THIS ASSUMES |S0| = 0 FOR ALL SUBSETS

clear all
clc

n = 20;
k = 4;

args = struct('n',n,'k',k,'type','kundir');

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
optx = intlinprog(c',[2:n+1],A,-b);
toc

disp('Max r-Robustness:')
optx(1)


% THIS DEFINITION IS OLD -- USE THE NEW DEFINITION WITH A STRUCT
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