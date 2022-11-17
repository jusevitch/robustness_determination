% Brute Force r-Robustness Verification
% Theory and Algorithm developed by James Usevitch

clear all
clc

% NOTE:
% This method is outdated in two ways:
% 1. Determining the max r-robustness can be formulated as a zero-one
% integer program, which is more efficient than this method
% 2. r-Robustness is NOT determined for all graphs by the case when
% S1 U S2 = V. This is proven in another one of my papers.
% Therefore, this program should not be used regularly.

% WARNING! Do not set n too high!
% !!! DO NOT USE THIS WITH LARGE NUMBERS OF n!!!
% Specifically, don't use with numbers much larger than n=22. Trying to
% create nchoosek(1:30,15) might crash a laptop.
n = 19; 
k = 18;

% Create the Laplacian of the network to be analyzed

% Complete graph
% D = Dmatrix(n,[],'arbcomplete',0);
% L = D*D';

% k-Circulant directed graph
% L = kCirculant(n,k,'dir');

% k-Circulant undirected graph
L = kCirculant(n,k,'undir');

% k-nearest neighbor platoons;
% L = zeros(n);
% for ii=1:1:k
%     L = L + diag(-ones(n - ii,1),ii) + diag(-ones(n-ii,1),-ii);
% end
% L = L + abs(diag(L*ones(n,1)));


nodes = 1:n;

% For integer programming implementation
% Mi = kron(eye(n),[1,i]);

% Initialize reachvar at -1, for coding purposes. The code seeks to find the smallest value r for
% which there exists an r-reachable node for the set tuple (S1,S2). This
% smallest value is the r-robustness, since every set will have a node with
% reachability greater than or equal to that value.

reachvar = -1;
lowestSvec = zeros(n,1);
lowestcertSet = [];

% Set a limit for how many lowest sets are saved
certSetlim = 20;


% Brute force computation of r-Robustness
% Don't use i for an indexing variable -- it needs to be the complex
% variable
for k=1:1:n-1
    if k ~= 1
        clear Sidx
    end
    
    Sidx = nchoosek(nodes,k);
    
    for jj =1:1:size(Sidx,1)
        
        Svec = ones(n,1);
        Svec(Sidx(jj,:)) = i;
        rhotemp = norm((1+i)/2*L*Svec, inf);
        
        % If rhotemp is less than reachvar, set reachvar equal to rhotemp.
        % Also return sets which has the lowest reachability
        % If not, reachvar remains the same.
        if reachvar == -1 || rhotemp < reachvar
            reachvar = rhotemp;
            lowestSvec = Svec;
            lowestcertSet = [Svec];
        elseif rhotemp == reachvar
            if size(lowestcertSet,2) < certSetlim
                lowestcertSet = [lowestcertSet Svec];
                % Change the sets saved
%             else
%                 lowestcertSet = [lowestcertSet]
            end
                    
        end
        
    end
    
end

disp('Calculated r-Robustness: ')
reachvar
