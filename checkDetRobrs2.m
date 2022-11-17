function outstruct = checkDetRobrs2(results)

% Checks difference between DetRob and rs2 

rdiffvec = find(round(results.DetRob(:,2)) - round(results.rs2(:,2)));

rposdiff = find(round(results.DetRob(:,2)) - round(results.rs2(:,2)) > 0);
rnegdiff = find(round(results.DetRob(:,2)) - round(results.rs2(:,2)) < 0);

sdiffvec = find(round(results.DetRob(:,5)) - round(results.rs2(:,5)));

sposdiff = find(round(results.DetRob(:,5)) - round(results.rs2(:,5)) > 0);
snegdiff = find(round(results.DetRob(:,5)) - round(results.rs2(:,5)) < 0);

diffvec = union(rdiffvec,sdiffvec); % Combines two vectors with no repeats

% Boolean vector that will tell whether the rsGnl2 results which are
% different than the DetRob results have a valid (S1,S2) certificate. 
certvec = zeros(length(diffvec),1); 

TM = results.testedMatrices(find(cell2mat(results.testedMatrices(:,end)) == 1),:);


for ii = 1:1:length(diffvec)
    
    TMii = TM(ii,:);
    
    L = TMii{1};
    
    tempstr = rsGnl2(struct('L',L));
    r = round(tempstr.r);
    s = round(tempstr.s);
    
    S1 = round(tempstr.two.b1);
    S2 = round(tempstr.two.b2);
    
    reachvec = [L*S1; L*S2];
    
    if length(reachvec(reachvec >= r)) >= s
        certvec(ii) = 1;
    end
    
end

outstruct.diffvec = diffvec;
outstruct.certvec = certvec;
outstruct.rposdiff = rposdiff;
outstruct.rnegdiff = rnegdiff;
outstruct.sposdiff = sposdiff;
outstruct.snegdiff = snegdiff;

end