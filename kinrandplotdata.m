function kinrandplotdata(results,fignumin)

% Plots the data in a results struct coming as output of rRobBenchmark2
% The results struct must have one or more of the following fields:
%   finalLB
%   finalintlin
%   finalintlinwm
%   finalgurobi
%   finalgurobiwm
%
% Optional fields:
%   results.plotalgs : vector with 1's for every algorithm's data to be plotted,
%                      and 0's otherwise. Algorithms are                    
%                      [LB intlin intlinwm gurobi gurobiwm]. If this field
%                      does not exist, all algorithms (with data in
%                      results) will be plotted.
%   results.separatefigs : set to 1 for all plots to have their own figures
%                          (no subplots used). Set to 0 otherwise. Default
%                          value is 0.
%
% Optional argument:
% fignumin : number of the figure to plot to.

unpack(results)

if isfield(results,'plotalgs')
    plotalgs = results.plotalgs;
else
    plotalgs = [1 1 1 1 1];
end

% Sets the figure number

if nargin == 2
    fignum = fignumin;
else
    fignum = 1;
end

% Sets a boolean depending on value of results.sepfigs

if isfield(results,'sepfigs') && results.sepfigs == 1
    sepfigs = 1;
else
    sepfigs = 0;
end

% Horizontal spacing of box plots

spc = .2;

% LB results

if plotalgs(1) && isfield(results,'finalLB')
    if isfield(finalLB,'kinrand')
        
        nvec = unique(finalLB.kinrand(:,3));
        kvec = unique(finalLB.kinrand(:,4));
        
        % Split data by p and n value; p is 1st index, n is 2nd index
        
        splitdata = cell(length(kvec),length(nvec));
        
        % Vectors for plotting computation time indexed by length of kvec
        
        meanvecs = zeros(length(kvec),length(nvec)); % rows relate to kvec, columns to nvec
        maxvecs = zeros(length(kvec),length(nvec));
        minvecs = zeros(length(kvec),length(nvec));
        
        for ii=1:1:length(kvec)
            for jj=1:1:length(nvec)
                k =  kvec(ii);
                n = nvec(jj);
                
                %         tempidx = find(finalLB.kinrand(:,4) == k & finalLB.kinrand(:,3) == n);
                
                splitdata{ii,jj} = finalLB.kinrand(finalLB.kinrand(:,4) == k & finalLB.kinrand(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(kvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else            
                subplot(1,length(kvec),kk);
            end
            
%             plot(nvec,meanvecs(kk,:)','^-b',nvec,maxvecs(kk,:),'--b',nvec,minvecs(kk,:),'-.b')
            errorbar(nvec,meanvecs(kk,:),(meanvecs(kk,:) - minvecs(kk,:)),(maxvecs(kk,:) - meanvecs(kk,:)),'o-b')
            
%             xlim([min(nvec) max(nvec)]);
            set(gca, 'YScale', 'log')
            title(['k-in Random Digraph with k =  ' ' ' num2str(kvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('LB max', 'LB mean', 'LB min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.kvec = kvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No kinrand results for finalLB')
    end
else
    disp('No finalLB results')
end


% intlin results

if plotalgs(2) && isfield(results,'finalintlin')
    if isfield(finalintlin,'kinrand')
        
        nvec = unique(finalintlin.kinrand(:,3));
        kvec = unique(finalintlin.kinrand(:,4));
        
        % Split data by p and n value; p is 1st index, n is 2nd index
        
        splitdata = cell(length(kvec),length(nvec));
        
        % Vectors for plotting computation time indexed by length of kvec
        
        meanvecs = zeros(length(kvec),length(nvec)); % rows relate to kvec, columns to nvec
        maxvecs = zeros(length(kvec),length(nvec));
        minvecs = zeros(length(kvec),length(nvec));
        
        for ii=1:1:length(kvec)
            for jj=1:1:length(nvec)
                k =  kvec(ii);
                n = nvec(jj);
                
                %         tempidx = find(finalintlin.kinrand(:,4) == k & finalintlin.kinrand(:,3) == n);
                
                splitdata{ii,jj} = finalintlin.kinrand(finalintlin.kinrand(:,4) == k & finalintlin.kinrand(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(kvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else            
                subplot(1,length(kvec),kk);
            end
            hold on
%             plot(nvec,meanvecs(kk,:)','v-r',nvec,maxvecs(kk,:),'--r',nvec,minvecs(kk,:),'-.r')
            errorbar(nvec + spc*ones(length(nvec),1),meanvecs(kk,:),(meanvecs(kk,:) - minvecs(kk,:)),(maxvecs(kk,:) - meanvecs(kk,:)),'o-r')
            hold on
            set(gca, 'YScale', 'log')
            title(['k-in Random Digraph with k =  ' ' ' num2str(kvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('intlin max', 'intlin mean', 'intlin min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.kvec = kvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No kinrand results for finalintlin')
    end
else
    disp('No finalintlin results')
end




% intlinwm results

if plotalgs(3) && isfield(results,'finalintlinwm')
    if isfield(finalintlinwm,'kinrand')
        
        nvec = unique(finalintlinwm.kinrand(:,3));
        kvec = unique(finalintlinwm.kinrand(:,4));
        
        % Split data by p and n value; p is 1st index, n is 2nd index
        
        splitdata = cell(length(kvec),length(nvec));
        
        % Vectors for plotting computation time indexed by length of kvec
        
        meanvecs = zeros(length(kvec),length(nvec)); % rows relate to kvec, columns to nvec
        maxvecs = zeros(length(kvec),length(nvec));
        minvecs = zeros(length(kvec),length(nvec));
        
        for ii=1:1:length(kvec)
            for jj=1:1:length(nvec)
                k =  kvec(ii);
                n = nvec(jj);
                
                %         tempidx = find(finalintlinwm.kinrand(:,4) == k & finalintlinwm.kinrand(:,3) == n);
                
                splitdata{ii,jj} = finalintlinwm.kinrand(finalintlinwm.kinrand(:,4) == k & finalintlinwm.kinrand(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(kvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else            
                subplot(1,length(kvec),kk);
            end
            hold on
            plot(nvec,meanvecs(kk,:)','v-m',nvec,maxvecs(kk,:),'--m',nvec,minvecs(kk,:),'-.m')
            hold on
            set(gca, 'YScale', 'log')
            title(['k-in Random Digraph with k =  ' ' ' num2str(kvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('intlinwm max', 'intlinwm mean', 'intlinwm min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.kvec = kvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No kinrand results for finalintlinwm')
    end
else
    disp('No finalintlinwm results')
end




% gurobi results

if plotalgs(4) && isfield(results,'finalgurobi')
    if isfield(finalgurobi,'kinrand')
        
        nvec = unique(finalgurobi.kinrand(:,3));
        kvec = unique(finalgurobi.kinrand(:,4));
        
        % Split data by p and n value; p is 1st index, n is 2nd index
        
        splitdata = cell(length(kvec),length(nvec));
        
        % Vectors for plotting computation time indexed by length of kvec
        
        meanvecs = zeros(length(kvec),length(nvec)); % rows relate to kvec, columns to nvec
        maxvecs = zeros(length(kvec),length(nvec));
        minvecs = zeros(length(kvec),length(nvec));
        
        for ii=1:1:length(kvec)
            for jj=1:1:length(nvec)
                k =  kvec(ii);
                n = nvec(jj);
                
                %         tempidx = find(finalgurobi.kinrand(:,4) == k & finalgurobi.kinrand(:,3) == n);
                
                splitdata{ii,jj} = finalgurobi.kinrand(finalgurobi.kinrand(:,4) == k & finalgurobi.kinrand(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(kvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else            
                subplot(1,length(kvec),kk);
            end
            hold on
            plot(nvec,meanvecs(kk,:)','v-g',nvec,maxvecs(kk,:),'--g',nvec,minvecs(kk,:),'-.g')
            hold on
            set(gca, 'YScale', 'log')
            title(['k-in Random Digraph with k =  ' ' ' num2str(kvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('gurobi max', 'gurobi mean', 'gurobi min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.kvec = kvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No kinrand results for finalgurobi')
    end
else
    disp('No finalgurobi results')
end



% gurobiwm results

if plotalgs(5) && isfield(results,'finalgurobiwm')
    if isfield(finalgurobiwm,'kinrand')
        
        nvec = unique(finalgurobiwm.kinrand(:,3));
        kvec = unique(finalgurobiwm.kinrand(:,4));
        
        % Split data by p and n value; p is 1st index, n is 2nd index
        
        splitdata = cell(length(kvec),length(nvec));
        
        % Vectors for plotting computation time indexed by length of kvec
        
        meanvecs = zeros(length(kvec),length(nvec)); % rows relate to kvec, columns to nvec
        maxvecs = zeros(length(kvec),length(nvec));
        minvecs = zeros(length(kvec),length(nvec));
        
        for ii=1:1:length(kvec)
            for jj=1:1:length(nvec)
                k =  kvec(ii);
                n = nvec(jj);
                
                %         tempidx = find(finalgurobiwm.kinrand(:,4) == k & finalgurobiwm.kinrand(:,3) == n);
                
                splitdata{ii,jj} = finalgurobiwm.kinrand(finalgurobiwm.kinrand(:,4) == k & finalgurobiwm.kinrand(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(kvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else            
                subplot(1,length(kvec),kk);
            end
            hold on
            plot(nvec,meanvecs(kk,:)','v-c',nvec,maxvecs(kk,:),'--c',nvec,minvecs(kk,:),'-.c')
            hold on
            set(gca, 'YScale', 'log')
            title(['k-in Random Digraph with k =  ' ' ' num2str(kvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('gurobiwm max', 'gurobiwm mean', 'gurobiwm min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.kvec = kvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No kinrand results for finalgurobiwm')
    end
else
    disp('No finalgurobiwm results')
end

end