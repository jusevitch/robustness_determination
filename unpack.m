function unpack(results)

% Unpacks the results struct of rRobBenchmark2 into the calling function
% To unpack into the base workspace, use 'base' instead of 'caller'

if isfield(results,'finalLB')
    assignin('caller','finalLB',results.finalLB);
end

if isfield(results,'finalintlin')
    assignin('caller','finalintlin',results.finalintlin);
end

if isfield(results,'finalintlinwm')
    assignin('caller','finalintlinwm',results.finalintlinwm);
end

if isfield(results,'finalgurobi')
    assignin('caller','finalgurobi',results.finalgurobi);
end

if isfield(results,'finalgurobiwm')
    assignin('caller','finalgurobiwm',results.finalgurobiwm);
end

if isfield(results,'testedMatrices')
    assignin('caller','testedMatrices',results.testedMatrices);
end

end