% r-Robustness exhaustive search method
% Theory and Code by James Usevitch

% WARNING: For a digraph of n nodes, this program will search through all
% possible ways to divide the graph into S_1, S_2, and S_0. THIS HAS
% EXPONENTIAL COMPLEXITY. RUN ON LARGE GRAPHS AT YOUR OWN RISK.

% Also, try this in C++. It should go a whole lot faster.

clear all;
clc;

% Parallel processing option - set to 1 to use a parfor loop
paropt = 1;

n = 13;
k = 4;

L = makegraph('kundir',n,k);

% Random matrix methodology
% Entries of Adjacency matrix chosen randomly
% A better method should be implemented that makes typical random graphs such as
% Barabasi-Albert, Erdos-Renyi, etc.
% Adj = randi([0,1],n,n);
% D = diag(Adj*ones(n,1));
% L = D - Adj;

Lmx = [L zeros(size(L)); zeros(size(L)) L];
nvec = (n-1)*ones(2*n,1);

sig1 = zeros(n,1);
sig2 = zeros(n,1);

% Initialize optimal values
maxrob = n; % Robustness can never exceed n/2, so any set obj value will be less than this
sig1opt = zeros(n,1);
sig2opt = zeros(n,1);

% For parallel processing: [maxrob, sig1opt, sig2opt];
sigparopt = zeros(1,2*n+1);
sigparopt(1,1) = n;



if paropt == 0
    
    disp('2^n-1 is equal to:')
    iternum = 2^n-1
    
    for s1=1:1:2^(n)-1
        sig1 = de2bi(s1,n);
        sig1 = sig1';
        
        z1 = find(sig1 == 0);
        n2 = length(z1);
        for s2=1:1:2^(n2)-1
            v2 = de2bi(s2,n2);
            sig2 = zeros(n,1);
            sig2(z1) = v2;
            if size(sig2,1) == 1
                sig2 = sig2';
            end
            
            % Calculate objective
            p = norm(Lmx*[sig1; sig2] + nvec,'inf') - (n-1);
            
            if p < maxrob
                sig1opt = sig1;
                sig2opt = sig2;
                maxrob = p;
            end
        end
        disp(strcat(num2str(s1),' / ',num2str(iternum)))
    end
    
    disp('Robustness:')
    maxrob
    
    disp('One possible minimal S1 / S2 pair:')
    disp('S1:')
    sig1opt'
    disp('S2:')
    sig2opt'
    
elseif paropt == 1
    
    sigparopt = zeros(1,2*n+1);
    sigparopt(1,1) = n;
    
    parfor s1=1:2^n-1
        
        sig1 = de2bi(s1,n)';
        %         sig1 = sig1';
        
        z1 = find(sig1 == 0);
        n2 = length(z1);
        for s2=1:1:2^(n2)-1
            v2 = de2bi(s2,n2);
            sig2 = zeros(n,1);
            sig2(z1) = v2;
            if size(sig2,1) == 1
                sig2 = sig2';
            end
            
            % Calculate objective
            p = norm(Lmx*[sig1; sig2] + nvec,'inf') - (n-1);
            
            sigparopt = compareSig(sigparopt, [p sig1' sig2']);
            %             if p < min(sigparopt(:,1))
            %                 %                 sig1opt = sig1;
            %                 %                 sig2opt = sig2;
            %                 sigparopt = [sigparopt; p sig1' sig2']
            %             end
        end
    end
    
    disp('Maximum r-Robustness:')
    sigparopt(1)
    
    sig1opt = sigparopt(1,2:n+1)';
    sig2opt = sigparopt(1,n+2:end)';
    
    disp('Sets: (One possible example)')
    disp('S1:')
    sig1opt'
    disp('S2:')
    sig2opt'
    
end




function outvec = compareSig(oldlist,newvec)

if isempty(oldlist)
    outvec = newvec;
elseif newvec(1,1) < oldlist(1,1)
    outvec = newvec;
else
    outvec = oldlist;
end

end


% The makegraph function
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