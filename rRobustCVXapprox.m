% r-Robust Linear Program Approximation
% Solves the r-Robust convex problem WITHOUT constraint that the 0-1 vector
% be integer.

clear all
clc

n = 50;
k = 3;

% Create the Laplacian of the network to be analyzed

L = makegraph('kplatoon',n,k);
Lk = makegraph('complete',n,k); % Complete graph Laplacian for constraints
Lk1 = Lk(1,:);

% Convex solution
% DOES NOT WORK. The solver simply converges towards a vector in span(ones) which is
% between 0 and 1. 

testvec = [1 zeros(1,n-1)]';

cvx_begin
    variable t nonnegative;
    variable x(n) nonnegative;
    variable x2(n) nonnegative;
    
    minimize t - sum_log((x - x2).^2)
    subject to
        -t*ones(n,1) <= L*x <= t*ones(n,1);
        x <= ones(n,1);
        x2 <= ones(n,1);
        x + x2 == ones(n,1);
cvx_end

t
x

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

% Old
% testvec = [1 zeros(1,n-1)]';
% 
% cvx_begin
%     variable t nonnegative;
%     variable x(n) nonnegative;
%     
%     minimize t
%     subject to
%         -t*ones(n,1) <= L*x <= t*ones(n,1);
%         x <= ones(n,1);
%         ones(1,n)*x >= 1;
%         (ones(n,1) - x)'*ones(n,1) >= 1;
%         testvec'*x == 1;
%     
% cvx_end