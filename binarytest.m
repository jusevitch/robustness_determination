
% This script demonstrates that when trying to create a binary
% representation of splitting a set into 2 subsets (1 = subset1, 0 =
% subset2), it is only necessary to count to 2^(n-1)-1.
%
% In other words, suppose you have n objects. The jth entry of a binary
% vector of length n represents whether object j is in set 1 (value of 1)
% or set 2 (value of 0). 
%
% There are 2^n-2 possible ways to set the vector such that neither set is
% empty. But after counting to 2^(n-1)-1, the sets repeat themselves. To
% see this, subtract set3 from set1 below. Or compare set1 and set2.

q = 4;

for ii=1:1:2^q-2
    set1(:,ii) = de2bifree(ii,q)';
    set2(:,ii) = ones(q,1) - set1(:,ii);    
end

set1
set2


for ii=1:1:2^(q-1)-1
    set3(:,ii) = de2bifree(ii,q)';
    set3(:,2^q-1-ii) = ones(q,1) - set3(:,ii);
    
end