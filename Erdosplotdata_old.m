function outstruct = Erdosplotdata_old(instruct)

% These are matrices, NOT structs
finalLB = instruct.finalLB.Erdos;
finalintlin = instruct.finalintlin.Erdos;
finalintlinwm = instruct.finalintlinwm.Erdos;

if isfield(instruct,'finalgurobi')
    testgurobi = 1;
    finalgurobi = instruct.finalgurobi.Erdos;
end

if isfield(instruct,'finalgurobiwm')
    testgurobi = 1;
    finalgurobiwm = instruct.finalgurobiwm.Erdos;
end

% Calculate nvec

nvec = unique(finalLB(:,3));

% Strip out zeros

zerovec = find(finalLB(:,2) == 0);
finalLB_z = finalLB(zerovec);
finalLB(zerovec,:) = [];

zerovec = find(finalintlin(:,2) == 0)
finalintlin_z = finalintlin(zerovec);
finalintlin(zerovec,:) = [];

zerovec = find(finalintlinwm(:,2) == 0)
finalintlinwm_z = finalintlinwm(zerovec);
finalintlinwm(zerovec,:) = [];

if(testgurobi)
    zerovec = find(finalgurobi(:,2) == 0)
    finalgurobi_z = finalgurobi(zerovec);
    finalgurobi(zerovec,:) = [];
    
    zerovec = find(finalgurobiwm(:,2) == 0)
    finalgurobiwm_z = finalgurobiwm(zerovec);
    finalgurobiwm(zerovec,:) = [];
end

outstruct.finalLB.max = 


% Calculate max, mean, min vectors



% Return




end