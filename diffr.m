function outstruct = diffr(results)

% Quickly compares the final r-robust results between algorithms tested
% with rRobBenchmark 2

% results1 is the first results structure, results2 is the second
% Returns a vector of indices of trials with differing max r-robust values

% Currently just tests Erdos, randdigraph, kinrand, koutrand



if isfield(results,'finalLB') && isfield(results,'finalintlin')
    
    graphset = {'Erdos','randdir','kinrand','koutrand'};
    
    for ii=1:1:length(graphset)
        
        if isfield(results.finalLB,graphset{ii}) && isfield(results.finalintlin,graphset{ii})
            
            if ~isempty(results.finalLB.(graphset{ii})) && ~isempty(results.finalLB.(graphset{ii}))
                
                outstruct.(graphset{ii}) = find(round(results.finalLB.(graphset{ii})(:,2)) - round(results.finalintlin.(graphset{ii})(:,2)));
                
                disp([graphset{ii} ' difference vector:'])
                outstruct.(graphset{ii})
                
            else
                error(['Either finalLB.' graphset{ii} ' or finalintlin.' graphset{ii} ' is empty.'])
                
            end
            
        end
        
    end
    
    %     if isfield(results.finalLB,'Erdos') && isfield(results.finalintlin,'Erdos')
    %
    %         if ~isempty(results.finalLB.Erdos) && ~isempty(results.finalLB.Erdos)
    %
    %             disp('Erdos difference vector:')
    %             outstruct.Erdos = find(round(results.finalLB.Erdos(:,2)) - round(results.finalintlin.Erdos(:,2)))
    %
    %         end
    %
    %     end
    
    
    
else
    error('Results not comparable--one or the other does not have either LB or intlin results')
end


end