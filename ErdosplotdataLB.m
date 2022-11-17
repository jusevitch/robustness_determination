function outstruct = ErdosplotdataLB(finalLB)

% Plots the Erdos portion of the finalLB data

nvec = unique(finalLB.Erdos(:,3));
pvec = unique(finalLB.Erdos(:,4));

% Split data by p and n value; p is 1st index, n is 2nd index

splitdata = cell(length(pvec),length(nvec));

% Vectors for plotting computation time indexed by length of pvec

meanvecs = zeros(length(pvec),length(nvec)); % rows relate to pvec, columns to nvec
maxvecs = zeros(length(pvec),length(nvec));
minvecs = zeros(length(pvec),length(nvec));

for ii=1:1:length(pvec)
    for jj=1:1:length(nvec)
        p = pvec(ii);
        n = nvec(jj);
        
%         tempidx = find(finalLB.Erdos(:,4) == p & finalLB.Erdos(:,3) == n);
        
        splitdata{ii,jj} = finalLB.Erdos(finalLB.Erdos(:,4) == p & finalLB.Erdos(:,3) == n,:);
        
        meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
        maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
        minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
        
    end
end

fig = figure(1)
for kk=1:1:length(pvec)
    
    subplot(1,length(pvec),kk)
    
    plot(nvec,meanvecs(kk,:)','^-b',nvec,maxvecs(kk,:),'--b',nvec,minvecs(kk,:),'-.b')

    set(gca, 'YScale', 'log')
    title(['p = ' ' ' num2str(pvec(kk))])
    xlabel('n values')
    ylabel('Max value of r')
    
end

outstruct.nvec = nvec;
outstruct.pvec = pvec;

outstruct.splitdata = splitdata;
outstruct.meanvecs = meanvecs;
outstruct.maxvecs = maxvecs;
outstruct.minvecs = minvecs;
outstruct.fig = fig;


end