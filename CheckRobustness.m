function outstruct = CheckRobustness(A,r,s)

% Implementation of CheckRobustness
% Theory by LeBlanc, Koutsoukos
% Code by James Usevitch

% REVISE THIS TO USE MORE EFFICIENT METHOD OF DetermineRobustness

n = size(A,1);
isRRobust = true;
outstruct.S1 = zeros(n,1);
outstruct.S2 = zeros(n,1);

for ii=1:1:2^(n)-1
    S1 = de2bifree(ii,n)';
    zervec = find(S1 == 0);
    for jj=1:1:2^(length(zervec))-1
        preS2 = de2bifree(jj,length(zervec));
        S2 = zeros(n,1);
        S2(zervec) = preS2;
        
        robust = RobustHolds(A,S1,S2,r,s);
        
        if robust == false % If S1, S2 found which prove graph is not (r,s)-robust for given values
            isRRobust = false;
            outstruct.isRRobust = isRRobust;
            outstruct.S1 = S1;
            outstruct.S2 = S2;
            return % Exits the function
        end
        
    end
    
end

% If the graph is (r,s)-robust, return the following values
outstruct.isRRobust = isRRobust;
outstruct.S1 = [];
outstruct.S2 = [];

end