function outstruct = rsCombine(incellmatrix)

% incellmatrix must be a cell matrix containing all results structs that
%   need to be combined.

names = fieldnames(incellmatrix{1});
outstruct = struct();

for ii=1:1:length(names)
    for jj=1:1:length(incellmatrix)
        if ~isfield(outstruct,names{ii})
            outstruct.(names{ii}) = incellmatrix{jj}.(names{ii});
        else
            outstruct.(names{ii}) = [outstruct.(names{ii}); incellmatrix{jj}.(names{ii})];
        end
    end
end


end