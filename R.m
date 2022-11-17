function outmulti = R(n)

% Returns R(n), the number of total possible unique pairs of nonempty, disjoint subsets
% in a graph with n agents

if isscalar(n)
    
    outmulti = 0;
    
    for kk = 2:1:n
        outmulti = outmulti + nchoosek(n,kk)*(2^(kk-1)-1);
    end
    
end

if ~isscalar(n) && isvector(n)
    outmulti = zeros(length(n),1);
    for ii=1:1:length(n)
        nn = n(ii);
        
        for kk = 2:1:nn
%                         outmulti(ii) = vpa(outmulti(ii) + nchoosek(nn,kk)*(2^(kk-1)-1),500);
            outmulti(ii) = outmulti(ii) + nchoosek(nn,kk)*(2^(kk-1)-1);
        end
    end
    
end

end