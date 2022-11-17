function randdirplotdata(results,fignumin)

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

if plotalgs(1) && isfield(results,'finalLB')
    if isfield(finalLB,'randdir')
        
        nvec = unique(finalLB.randdir(:,3));
        pvec = unique(finalLB.randdir(:,4));
        
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
                
                %         tempidx = find(finalLB.randdir(:,4) == p & finalLB.randdir(:,3) == n);
                
                splitdata{ii,jj} = finalLB.randdir(finalLB.randdir(:,4) == p & finalLB.randdir(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(pvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else
                subplot(1,length(pvec),kk);
            end
            
            
            
            %             plot(nvec,meanvecs(kk,:)','^-b',nvec,maxvecs(kk,:),'--b',nvec,minvecs(kk,:),'-.b')
            errorbar(nvec,meanvecs(kk,:),(meanvecs(kk,:) - minvecs(kk,:)),(maxvecs(kk,:) - meanvecs(kk,:)),'o-b')
            
            set(gca, 'YScale', 'log')
            title(['Random Digraph with p = ' ' ' num2str(pvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('LB max', 'LB mean', 'LB min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.pvec = pvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No randdir results for finalLB')
    end
else
    disp('No finalLB results')
end


% intlin portion

if plotalgs(2) && isfield(results,'finalintlin')
    if isfield(finalintlin,'randdir')
        
        nvec = unique(finalintlin.randdir(:,3));
        pvec = unique(finalintlin.randdir(:,4));
        
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
                
                %         tempidx = find(finalintlin.randdir(:,4) == p & finalintlin.randdir(:,3) == n);
                
                splitdata{ii,jj} = finalintlin.randdir(finalintlin.randdir(:,4) == p & finalintlin.randdir(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(pvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else
                subplot(1,length(pvec),kk);
            end
            hold on
%             plot(nvec,meanvecs(kk,:)','v-r',nvec,maxvecs(kk,:),'--r',nvec,minvecs(kk,:),'-.r')
            errorbar(nvec + spc*ones(length(nvec),1),meanvecs(kk,:),(meanvecs(kk,:) - minvecs(kk,:)),(maxvecs(kk,:) - meanvecs(kk,:)),'o-r')
            hold off
            set(gca, 'YScale', 'log')
            title(['Random Digraph with p = ' ' ' num2str(pvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('LB max', 'LB mean', 'LB min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.pvec = pvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No randdir results for finalintlin')
    end
else
    disp('No finalintlin results')
end



% intlinwm portion

if plotalgs(3) && isfield(results,'finalintlinwm')
    if isfield(finalintlinwm,'randdir')
        
        nvec = unique(finalintlinwm.randdir(:,3));
        pvec = unique(finalintlinwm.randdir(:,4));
        
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
                
                %         tempidx = find(finalintlinwm.randdir(:,4) == p & finalintlinwm.randdir(:,3) == n);
                
                splitdata{ii,jj} = finalintlinwm.randdir(finalintlinwm.randdir(:,4) == p & finalintlinwm.randdir(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(pvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else
                subplot(1,length(pvec),kk);
            end
            hold on
            plot(nvec,meanvecs(kk,:)','v-m',nvec,maxvecs(kk,:),'--m',nvec,minvecs(kk,:),'-.m')
            hold off
            set(gca, 'YScale', 'log')
            title(['Random Digraph with p = ' ' ' num2str(pvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('LB max', 'LB mean', 'LB min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.pvec = pvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No randdir results for finalintlinwm')
    end
else
    disp('No finalintlinwm results')
end



% gurobi portion

if plotalgs(4) && isfield(results,'finalgurobi')
    if isfield(finalgurobi,'randdir')
        
        nvec = unique(finalgurobi.randdir(:,3));
        pvec = unique(finalgurobi.randdir(:,4));
        
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
                
                %         tempidx = find(finalgurobi.randdir(:,4) == p & finalgurobi.randdir(:,3) == n);
                
                splitdata{ii,jj} = finalgurobi.randdir(finalgurobi.randdir(:,4) == p & finalgurobi.randdir(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(pvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else
                subplot(1,length(pvec),kk);
            end
            hold on
            plot(nvec,meanvecs(kk,:)','v-g',nvec,maxvecs(kk,:),'--g',nvec,minvecs(kk,:),'-.g')
            hold off
            set(gca, 'YScale', 'log')
            title(['Random Digraph with p = ' ' ' num2str(pvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('LB max', 'LB mean', 'LB min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.pvec = pvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No randdir results for finalgurobi')
    end
else
    disp('No finalgurobi results')
end



% gurobiwm portion

if plotalgs(5) && isfield(results,'finalgurobiwm')
    if isfield(finalgurobiwm,'randdir')
        
        nvec = unique(finalgurobiwm.randdir(:,3));
        pvec = unique(finalgurobiwm.randdir(:,4));
        
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
                
                %         tempidx = find(finalgurobiwm.randdir(:,4) == p & finalgurobiwm.randdir(:,3) == n);
                
                splitdata{ii,jj} = finalgurobiwm.randdir(finalgurobiwm.randdir(:,4) == p & finalgurobiwm.randdir(:,3) == n,:);
                
                meanvecs(ii,jj) = mean(splitdata{ii,jj}(:,1));
                maxvecs(ii,jj) = max(splitdata{ii,jj}(:,1));
                minvecs(ii,jj) = min(splitdata{ii,jj}(:,1));
                
            end
        end
        
        if ~sepfigs
            fig = figure(fignum); % Single figure if sepfigs = 0
        end
        
        for kk=1:1:length(pvec)
            
            if sepfigs
                fig{kk} = figure(fignum + (kk-1)); % Plot in separate figures (fignum, fignum+1, fignum+2)
            else
                subplot(1,length(pvec),kk);
            end
            hold on
            plot(nvec,meanvecs(kk,:)','v-c',nvec,maxvecs(kk,:),'--c',nvec,minvecs(kk,:),'-.c')
            hold off
            set(gca, 'YScale', 'log')
            title(['Random Digraph with p = ' ' ' num2str(pvec(kk))])
            xlabel('n values')
            ylabel('Computation Time (sec)')
            %     legend('LB max', 'LB mean', 'LB min')
            
        end
        
        outstruct.nvec = nvec;
        outstruct.pvec = pvec;
        
        outstruct.splitdata = splitdata;
        outstruct.meanvecs = meanvecs;
        outstruct.maxvecs = maxvecs;
        outstruct.minvecs = minvecs;
        outstruct.fig = fig;
        
    else
        disp('No randdir results for finalgurobiwm')
    end
else
    disp('No finalgurobiwm results')
end

end