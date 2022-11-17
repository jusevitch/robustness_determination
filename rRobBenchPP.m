% Postprocessing for rRobBenchmark

% This file plots the data from the rRobBenchmark.m file

nvec = args.nvec;
pvec = args.pvec;

% Erdos-Renyi results

LBcell = {};
intlincell = {};
intlinwmcell = {};

if testGurobi
    gurobicell = {};
    gurobiwmcell = {};
end

for ii=1:1:length(pvec)
    
    idx = find(finalLB.Erdos(:,4) == pvec(ii));
    LBcell{ii} = finalLB.Erdos(idx,:);
    
    idx = find(finalintlin.Erdos(:,4) == pvec(ii));
    intlincell{ii} = finalintlin.Erdos(idx,:);
    
    idx = find(finalintlinwm.Erdos(:,4) == pvec(ii));
    intlinwmcell{ii} = finalintlinwm.Erdos(idx,:);
    
    if testGurobi
        idx = find(finalgurobi.Erdos(:,4) == pvec(ii));
        gurobicell{ii} = finalgurobi.Erdos(idx,:);
        
        idx = find(finalgurobiwm.Erdos(:,4) == pvec(ii));
        gurobiwmcell{ii} = finalgurobiwm.Erdos(idx,:);
    end
end

% r-Robustness values
fig1 = figure(1);
fig1.Name = 'Erdos-Renyi - Calculated r-Robustness Values';
clf
for ii=1:1:length(pvec)
    
    subplot(1,3,ii)
    title(['p value:' ' ' num2str(pvec(ii))]);
    hold on
    
    gLB = plot(nvec,round(LBcell{ii}(:,2)),'^');
    gi = plot(nvec,round(intlincell{ii}(:,2)),'v');
    giwm = plot(nvec,round(intlinwmcell{ii}(:,2)),'x');
    
    if testGurobi
        gg = plot(nvec,round(gurobicell{ii}(:,2)),'o');
        ggwm = plot(nvec,round(gurobiwmcell{ii}(:,2)),'s');
    end
    
    xlabel('Value of n')
    ylabel('r-Robustness')
    hold off
end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end

% Computation time

fig2 = figure(2);
fig2.Name = 'Erdos-Renyi - Computation Time';
clf
for ii=1:1:length(pvec)
    
    subplot(1,3,ii)
    title(['p value:' ' ' num2str(pvec(ii))]);
    %     ylim([0 .3])
    
    hold on
    gLB = plot(nvec,LBcell{ii}(:,1),'-^');
    gi = plot(nvec,intlincell{ii}(:,1),'-v');
    giwm = plot(nvec,intlinwmcell{ii}(:,1),'-x');
    
    if testGurobi
        gg = plot(nvec,gurobicell{ii}(:,1),'-o');
        ggwm = plot(nvec,gurobiwmcell{ii}(:,1),'-s');
    end
    
    xlabel('Value of n')
    ylabel('Computation time (sec)')
    hold off
    set(gca, 'YScale', 'log')
end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end


% Complete Graph Results

% r-Robustness values
fig3 = figure(3);
fig3.Name = 'Complete Graphs - Calculated r-Robustness Values';
clf
% for ii=1:1:length(pvec)

%     subplot(1,3,ii)
%     title(['p value:' ' ' num2str(pvec(ii))]);
hold on

gLB = plot(nvec,round(finalLB.complete(:,2)),'^');
gi = plot(nvec,round(finalintlin.complete(:,2)),'v');
giwm = plot(nvec,round(finalintlinwm.complete(:,2)),'x');

if testGurobi
    gg = plot(nvec,round(finalgurobi.complete(:,2)),'o');
    ggwm = plot(nvec,round(finalgurobiwm.complete(:,2)),'s');
end

xlabel('Value of n')
ylabel('r-Robustness')
hold off
% end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end

% Computation time

fig4 = figure(4);
fig4.Name = 'Complete Graphs - Computation Time';
clf
% for ii=1:1:length(pvec)

%     subplot(1,3,ii)
%     title(['p value:' ' ' num2str(pvec(ii))]);
%     ylim([0 .3])

hold on
gLB = plot(nvec,finalLB.complete(:,1),'^-');
gi = plot(nvec,finalintlin.complete(:,1),'v-');
giwm = plot(nvec,finalintlinwm.complete(:,1),'x-');

if testGurobi
    gg = plot(nvec,finalgurobi.complete(:,1),'o-');
    ggwm = plot(nvec,finalgurobiwm.complete(:,1),'s-');
end

xlabel('Value of n')
ylabel('Computation time (sec)')
hold off
set(gca, 'YScale', 'log')
% end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end


% k-Circulant Undirected

LBcell = {};
intlincell = {};
intlinwmcell = {};

if testGurobi
    gurobicell = {};
    gurobiwmcell = {};
end


kvec = 2:3:floor(20/2); % which is max value of n-2

% ATTENTION: The following code addresses a goof I made with the
% rRobBenchmark. Fix it later.

for ii=1:1:length(kvec)
    
    idx = find(finalLB.kundir(:,4) == kvec(ii));
    LBcell{ii} = finalLB.kundir(idx,:);
    
    idx = find(finalintlin.kundir(:,4) == kvec(ii));
    intlincell{ii} = finalintlin.kundir(idx,:);
    
    idx = find(finalintlinwm.kundir(:,4) == kvec(ii));
    intlinwmcell{ii} = finalintlinwm.kundir(idx,:);
    
    if testGurobi
        idx = find(finalgurobi.kundir(:,4) == kvec(ii));
        gurobicell{ii} = finalgurobi.kundir(idx,:);
        
        idx = find(finalgurobiwm.kundir(:,4) == kvec(ii));
        gurobiwmcell{ii} = finalgurobiwm.kundir(idx,:);
    end
end

% r-Robustness values: n vs k
fig5 = figure(5);
fig5.Name = 'k-Circulant Undirected - Calculated r-Robustness Values - varying n';
clf
for ii=1:1:length(kvec)
    
    subplot(1,3,ii)
    title(['k value:' ' ' num2str(kvec(ii))]);
    hold on
    
    gLB = plot(LBcell{ii}(:,3),round(LBcell{ii}(:,2)),'^');
    gi = plot(intlincell{ii}(:,3),round(intlincell{ii}(:,2)),'v');
    giwm = plot(intlinwmcell{ii}(:,3),round(intlinwmcell{ii}(:,2)),'x');
    
    if testGurobi
        gg = plot(gurobicell{ii}(:,3),round(gurobicell{ii}(:,2)),'o');
        ggwm = plot(gurobiwmcell{ii}(:,3),round(gurobiwmcell{ii}(:,2)),'s');
    end
    
    
    xlabel('Value of n')
    ylabel('r-Robustness')
    hold off
    ylim([0 10])
end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end

% Computation Time n vs k

fig6 = figure(6);
fig6.Name = 'k-Circulant Undirected - Computation Time - varying n';
clf
for ii=1:1:length(kvec)
    
    subplot(1,3,ii)
    title(['k value:' ' ' num2str(kvec(ii))]);
    hold on
    
    gLB = plot(LBcell{ii}(:,3),LBcell{ii}(:,1),'-^');
    gi = plot(intlincell{ii}(:,3),intlincell{ii}(:,1),'-v');
    giwm = plot(intlinwmcell{ii}(:,3),intlinwmcell{ii}(:,1),'-x');
    
    if testGurobi
        gg = plot(gurobicell{ii}(:,3),gurobicell{ii}(:,1),'-o');
        ggwm = plot(gurobiwmcell{ii}(:,3),gurobiwmcell{ii}(:,1),'-s');
    end
    
    xlabel('Value of n')
    ylabel('Computation time (sec)')
    hold off
    set(gca, 'YScale', 'log') % Sets y axis to be log scale
    ylim([1e-3 1e5])
end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end

% r-robustness values: k vs n

nplotvec = [13 17 20];

LBcell = {};
intlincell = {};
intlinwmcell = {};

if testGurobi
    gurobicell = {};
    gurobiwmcell = {};
end

for ii=1:1:length(nplotvec)
    
    idx = find(finalLB.kundir(:,3) == nplotvec(ii));
    LBcell{ii} = finalLB.kundir(idx,:);
    
    idx = find(finalintlin.kundir(:,3) == nplotvec(ii));
    intlincell{ii} = finalintlin.kundir(idx,:);
    
    idx = find(finalintlinwm.kundir(:,3) == nplotvec(ii));
    intlinwmcell{ii} = finalintlinwm.kundir(idx,:);
    
    if testGurobi
        idx = find(finalgurobi.kundir(:,3) == nplotvec(ii));
        gurobicell{ii} = finalgurobi.kundir(idx,:);
        
        idx = find(finalgurobiwm.kundir(:,3) == nplotvec(ii));
        gurobiwmcell{ii} = finalgurobiwm.kundir(idx,:);
    end
end

fig7 = figure(7);
fig7.Name = 'k-Circulant Undirected - Calculated r-Robustness Values - varying k';
clf
for ii=1:1:length(nplotvec)
    
    subplot(1,3,ii)
    title(['n value:' ' ' num2str(nplotvec(ii))]);
    hold on
    
    gLB = plot(kvec(1:size(LBcell{ii},1)),round(LBcell{ii}(:,2)),'^');
    gi = plot(kvec(1:size(intlincell{ii},1)),round(intlincell{ii}(:,2)),'v');
    giwm = plot(kvec(1:size(intlinwmcell{ii},1)),round(intlinwmcell{ii}(:,2)),'x');
    
    if testGurobi
        gg = plot(kvec(1:size(gurobicell{ii},1)),round(gurobicell{ii}(:,2)),'o');
        ggwm = plot(kvec(1:size(gurobiwmcell{ii},1)),round(gurobiwmcell{ii}(:,2)),'s');
    end
    
%     ylim = 
    xlabel('Value of n')
    ylabel('r-Robustness')
    hold off
end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end

fig8 = figure(8);
fig8.Name = 'k-Circulant Undirected - Computation Time - varying n';
clf
for ii=1:1:length(kvec)
    
    subplot(1,3,ii)
    title(['n value:' ' ' num2str(nplotvec(ii))]);
    hold on
    
    gLB = plot(kvec(1:size(LBcell{ii},1)),LBcell{ii}(:,1),'^-');
    gi = plot(kvec(1:size(intlincell{ii},1)),intlincell{ii}(:,1),'v-');
    giwm = plot(kvec(1:size(intlinwmcell{ii},1)),intlinwmcell{ii}(:,1),'x-');
    
    if testGurobi
        gg = plot(kvec(1:size(gurobicell{ii},1)),gurobicell{ii}(:,1),'o-');
        ggwm = plot(kvec(1:size(gurobiwmcell{ii},1)),gurobiwmcell{ii}(:,1),'s-');
    end
    
    xlabel('Value of k')
    ylabel('Computation time (sec)')
    hold off
    set(gca, 'YScale', 'log') % Sets y axis to be log scale
end

if testGurobi
    legend([gLB, gi, giwm, gg, ggwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start','Gurobi on 12 cores', 'Gurobi on 12 cores with warm start');
else
    legend([gLB, gi, giwm],'Nominal Algorithm','inlinprog on 1 core','intlinprog with warm start');
end

