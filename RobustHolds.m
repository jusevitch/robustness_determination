%#coder

% Implementation of RobustHolds
% Note: this code has slightly vectorized and altered the original proposed
% algorithm by LeBlanc to speed up implementation
% Code by James Usevitch

function isRRobust = RobustHolds(A,S1,S2,r,s)

% Input: Each Sk should be a binary vector of length n where Sk(j) = 1 if
% j \in Sk
% varargin must be {r} or {r,s}.
% n = size(A,1);

isRRobust = false;
% sr = [0; 0]; % sr1 and sr2
% bigS = {S1; S2}; % Cell matrix with S1 and S2
% r = varargin{1};
% if length(varargin) == 2
%     s = varargin{2};
% elseif length(varargin) == 1
%     s = 1; % r robustness is (r,1)-robustness
% end

% for k=1:1:2
%     reachvec = A(find(bigS{k}),:)*(ones(n,1) - bigS{k}); % Reachability of nodes in Sk; if Sk has 1's for nodes in Sk,
                                                         % (ones - Sk) will have 1's for nodes NOT in Sk.

reachvec = sum(A(S1 == 1,S1 == 0),2); % Reachability of nodes in Sk; if Sk has 1's for nodes in Sk, the vector
                                                         % (Sk == 0) will have 1's for nodes NOT in Sk.
sr1 = sum(reachvec >= r); % Number of nodes with reachability greater than r
if (sr1 == sum(S1))
    isRRobust = true;
    return
end
% end
reachvec = sum(A(S2 == 1,S2 == 0),2);% Reachability of nodes in Sk; if Sk has 1's for nodes in Sk,
                                                         % (ones - Sk) will have 1's for nodes NOT in Sk.
sr2 = sum(reachvec >= r); % Number of nodes with reachability greater than r
if (sr2 == sum(S2))
    isRRobust = true;
    return
end

% if (sr(1) == length(find(S1))) || (sr(2) == length(find(S2))) || (sr(1) + sr(2) >= s)
% if (sr1 == sum(S1)) || (sr2 == sum(S2)) || (sr1 + sr2 >= s)
if (sr1 + sr2 >= s)
    isRRobust = true;
    return
end

end