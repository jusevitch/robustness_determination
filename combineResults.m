function outstruct = combineResults(incell)

% Combines the result structs from several runs of rRobBenchmark2
%
%

alglist = {'finalLB','finalintlin','finalintlinwm','finalgurobi','finalgurobiwm'};

graphlist = {'Erdos','complete','kundir','kdir','randdir','kinrand','koutrand'};


for ii=1:1:length(alglist)
    for jj=1:1:length(graphlist)
        
        tempmx = [];
        
        for kk=1:1:length(incell)
            if isfield(incell{kk},alglist{ii}) && isfield(incell{kk}.(alglist{ii}),graphlist{jj})
                tempmx = [tempmx; incell{kk}.(alglist{ii}).(graphlist{jj})];
            end
        end
        
        if ~isempty(tempmx)
            tempmx = sortrows(tempmx,[3 2 1]);
        end
        
        outstruct.(alglist{ii}).(graphlist{jj}) = tempmx;
        
    end
end


end