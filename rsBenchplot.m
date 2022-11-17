function rsBenchplot(args,instruct)

% THIS IS ONLY COMPATIBLE WITH rsBench3 DATA
% FOR DATA CREATED BY rsBench4, USE rsBenchplot2.m

% args should be the args vector used to create the data and should
% contain:
%   args.nvec
%   args.pvec
%   args.kvec

% algorithms:
%   1: LB
%   2: DetRob
%   3: rRG
%   4: rRG2
%   5: rs2
%   6: Fmax
%   7: rlow
%   8: rup

algnames = {'LB','DetRob','rRG','rRG2','rs2','Fmax','rlow','rup'};
algfullnames = {'Mod. Det. Rob.','Det. Rob.','asdf','r-Robust MILP','(r,s)-Robust MILP','Fmax MILP','r-Rob. Lower Bnd','r-Rob. Upper Bnd'};
graphnames = {'Erdos-Renyi', 'Random Digraph', 'k-In Random Digraph', 'k-Out Random Digraph'};

if isfield(instruct,'LB')
    algs{1} = instruct.LB;
end
if isfield(instruct,'DetRob')
    algs{2} = instruct.DetRob;
end
if isfield(instruct,'rRG')
    algs{3} = instruct.rRG;
end
if isfield(instruct,'rRG2')
    algs{4} = instruct.rRG2;
end
if isfield(instruct,'rs2')
    algs{5} = instruct.rs2;
end
if isfield(instruct,'Fmax')
    algs{6} = instruct.Fmax;
end
if isfield(instruct,'rlow')
    algs{7} = instruct.rlow;
end
if isfield(instruct,'rup')
    algs{8} = instruct.rup;
end

% LB = instruct.LB;
% DetRob = instruct.DetRob;
% rRG = instruct.rRG;
% rRG2 = instruct.rRG2;
% rs2 = instruct.rs2;
% Fmax = instruct.Fmax;
% rlow = instruct.rlow;
% rup = instruct.rup;

pvec = args.pvec;
kvec = args.kvec;

fignum = 1;

plotcolors = distinguishable_colors(length(algs));


% Plot DetRob vs rs2

% 1: Erdos
% 2: randdir
% 3: kinrand
% 4: koutrand

close all

% DetermineRobustness and rs2
% Iterate through Graphs first

algplotlist = [];

for ii=1:4 % Graph
    
    if ii <= 2
        pkvec = pvec;
    else
        pkvec = kvec;
    end
    
    for kk=1:1:length(pkvec) % pkvec
        
        for jj=[2 5] % Algorithms (DetRob)
            
            if ~isempty(algs{jj})
                algplotlist = [algplotlist jj];
                data = algs{jj}(algs{jj}(:,end) == ii,:);
                
                pkdata = data(data(:,4) == pkvec(kk),:);
                
                nvec = unique(pkdata(:,3));
                
                for ll=1:length(nvec)
                    ntime_mat = pkdata(pkdata(:,3) == nvec(ll),:);
                    avgtime(ll) = mean(ntime_mat(:,1));
                    maxtime(ll) = max(ntime_mat(:,1));
                    mintime(ll) = min(ntime_mat(:,1));
                end
                
                figure(fignum)
                hold on
                errorbar(nvec,avgtime,(avgtime-mintime),(maxtime-avgtime),'o-','color',plotcolors(jj,:))
                hold off
                set(gca, 'YScale', 'log')
                xlabel('Number of nodes')
                ylabel('Computation time')
                title(['Graph for ' graphnames(ii) ' graph and pk=' num2str(pkvec(kk))])
                
                %                 fignum = fignum +1;
            end
        end
        legend(algfullnames(algplotlist));
        fignum = fignum + 1;
    end
    
    %     fignum = (length(pkvec)-1) + 1;
end


% LB, rRobustGnl2, rRobupper, rRoblower
algplotlist = [];

for ii=1:4
    
    if ii <= 2
        pkvec = pvec;
    else
        pkvec = kvec;
    end
    
    for kk=1:1:length(pkvec) % pkvec
        
        for jj=[1 4 7 8] % Algorithms
            
            if ~isempty(algs{jj})
                algplotlist = [algplotlist jj];
                data = algs{jj}(algs{jj}(:,end) == ii,:);
                
                pkdata = data(data(:,4) == pkvec(kk),:);
                
                nvec = unique(pkdata(:,3));
                
                for ll=1:length(nvec)
                    ntime_mat = pkdata(pkdata(:,3) == nvec(ll),:);
                    avgtime(ll) = mean(ntime_mat(:,1));
                    maxtime(ll) = max(ntime_mat(:,1));
                    mintime(ll) = min(ntime_mat(:,1));
                end
                
                figure(fignum)
                hold on
                errorbar(nvec,avgtime,(avgtime-mintime),(maxtime-avgtime),'o-','color',plotcolors(jj,:))
                hold off
                set(gca, 'YScale', 'log')
                xlabel('Number of nodes')
                ylabel('Computation time')
                title(['Graph for ' graphnames(ii) ' graph and pk=' num2str(pkvec(kk))])
                
                %                 fignum = fignum +1;
            end
        end
        
        legend(algfullnames(algplotlist));
        fignum = fignum + 1;
    end
    
    
end


% Plot the differences between objective values of rRobustGnl2, rRoblower,
% and rRob upper





end