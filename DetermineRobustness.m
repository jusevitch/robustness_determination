%#coder

function outstruct = DetermineRobustness(args)

% Implementation of DetermineRobustness
% Theory by Leblanc, Koutsoukos
% Code by James Usevitch
% Arguments:
%   args.A : adjacency matrix of graph
%   args.smax : integer upper bound of values of s for which you want to
%               check. For example, to check only for r-robustness (which is equivalent
%               to (r,1)-robustness), set args.smax = 1
%   args.q : Testing

if isfield(args,'A')
    A = args.A;
elseif isfield(args,'L')
    A = diag(diag(args.L)) - args.L;
else
    error('No adjacency or Laplacian matrix specified in args')
end

n = size(A,1);
indegvec = A*ones(n,1); % Vector of in-degrees of all nodes
% outdegvec = ones(1,n)*A; % Vector of out-degrees of all nodes
minindeg = min(indegvec); % Minimum in-degree of all nodes
% minoutdeg = min(outdegvec);

time = 0;
tic;

% Sets an upper bound on the value of r.
% If minimum in-degree of nodes in the graph is strictly positive, then
% the upper bound on r is the lesser of the min in-degree and ceil(n/2).
% If minimum in-degree of nodes in the graph is 0, then the upper bound on
% r is 1 (e.g. consider a rooted outbranching; rooted
% outbranching exists iff graph is 1-robust)
r = min(max(minindeg,1),ceil(n/2));

% Set an upper bound on values of s to consider, e.g. only testing for
% (r,1)-robustness (equivalent to r-robustness)
if isfield(args,'smax')
    s = args.smax;
else
    s = n;
end

setsChecked = 0;

% Total # of unique, nonempty, disjoint subset pairs to check
totalsets = 0;
for kk = 2:1:n
    totalsets = totalsets + nchoosek(n,kk)*(2^(kk-1)-1);
end

% There are two binary iterators in the next section
% The first binary iterator encodes all possible ways to choose a subset
% from the graph, effectivey creating two "meta sets" M0 and M1.
%
% The second binary iterator encodes all ways to split M1 into two
% nonempty, disjoint subsets S1 and S2. S0 is effectively represented by
% M0.
%
% This is more efficient as it does not repeat combinations of S1 and S2.
% See  See "binarytest.m" for an example.

powvec = 0:-1:1-n;

for ii=1:1:2^(n)-1
%     M1 = de2bifree(ii,n)'; % Choose the subset M1
%     M1 = rem(floor(ii*pow2(0:-1:1-n)),2)'; % Does what de2bifree does, but faster since not external function and no 'char' call.
    M1 = rem(floor(ii*pow2(powvec)),2)';
%     M1 = mod(floor(ii*pow2(0:-1:1-n)),2) % Is this faster?
    if sum(M1) >= 2 % Rules out the cases when |M1| = 0, since we want n choose k, k >= 2
        M1vec = (M1 == 1); % Makes M1vec a logical binary vector
        sumM1 = sum(M1vec);
        for jj=1:1:2^(sumM1-1)-1
%             preS1 = de2bifree(jj,length(M1vec))';
            preS1 = rem(floor(jj*pow2(0:-1:1-sumM1)),2)';
%             preS2 = ones(sumM1,1) - preS1;
            preS2 = (preS1 == 0);
            S1 = zeros(n,1);
            S2 = zeros(n,1);
            S1(M1vec) = preS1;
            S2(M1vec) = preS2;
            
            % testing
%             if (norm(S1 - [1;1;0;1;0;1;0]) == 0 && norm(S2 - [0;0;1;0;1;0;1]) == 0) || (norm(S1 - [0;0;1;0;1;0;1]) == 0 && norm(S2 - [1;1;0;1;0;1;0]) == 0)
%                 disp('stop')
%             end
            % end testing
            
            isRSRobust = RobustHolds(A,S1,S2,r,s);
            
            if (isRSRobust == false) && (s > 0) % If S1, S2 found which prove graph is not (r,s)-robust for given values
                s = s - 1;
            end
            %testing
%             if r == 1 && s == 1
%                 disp('stop')
%             end
            %end testing
            while (isRSRobust == false) && (r > 0)
                while (isRSRobust == false) && (s > 0)
                    isRSRobust = RobustHolds(A,S1,S2,r,s);
                    % testing
%                     if r == 1 && s == 1 %&& isRSRobust == false
%                         disp('stop')
%                     end
                    %end testing
                    if isRSRobust == false
                        s = s-1;
                    end
                end
                if isRSRobust == false
                    r = r-1;
                    if isfield(args,'smax')
                        s = args.smax;
                    else
                        s = n;
                    end
                end
            end
            
            if mod(setsChecked,1000000) == 0
                disp([num2str(setsChecked) ' / ' num2str(totalsets) ' (' num2str(setsChecked/totalsets*100) '% complete)'])
                if isfield(args,'q')
                    tempstring = [num2str(setsChecked) ' / ' num2str(totalsets) ' (' num2str(setsChecked/totalsets*100) '% complete)'];
                    send(args.q, sprintf('%s',tempstring));
                end
            end
            
            if r == 0
                time = toc;
                
                outstruct.r = r;
                outstruct.s = 0;
                outstruct.setsChecked = setsChecked;
                outstruct.time = time;
                
                return
            end
            
            setsChecked = setsChecked + 1;
        end
    end
end
time = toc;

outstruct.r = r;
outstruct.s = s;
outstruct.setsChecked = setsChecked;
outstruct.time = time;



end
