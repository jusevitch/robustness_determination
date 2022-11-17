function rsBenchplot2(instruct)

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
algfullnames = {'Mod. Det. Rob.','Det. Rob.','asdf','r-Rob. MILP','(r,s)-Rob. MILP','Fmax MILP','r-Rob. Lower Bnd','r-Rob. Upper Bnd'};
graphnames = {'Erdos-Renyi Graphs', 'Random Digraphs', 'k-In Random Digraphs', 'k-Out Random Digraphs'};

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

fignum = 1;

plotcolors = distinguishable_colors(length(algs));
plotcolors = [plotcolors(1,:); plotcolors(3,:); plotcolors(2,:); plotcolors(4:end,:)];
plotcolors = [plotcolors(3,:); plotcolors(1,:); plotcolors(2,:); plotcolors(4:end,:)];


% Plot DetRob vs rs2

% 1: Erdos
% 2: randdir
% 3: kinrand
% 4: koutrand

total_pkvec = cell(2,1);

for ii=1:4
    for mm=1:length(algs)
        if ~isempty(algs{mm})
            if ii == 1 || ii == 2
                total_pkvec{1} = union(total_pkvec{1},algs{mm}(algs{mm}(:,6) == ii,5));
            elseif ii ==3 || ii == 4
                total_pkvec{2} = union(total_pkvec{2},algs{mm}(algs{mm}(:,6) == ii,5));
            end
        end
    end
end

close all

% DetermineRobustness and rs2
% Iterate through Graphs first

algplotlist = [];

for ii=1:4 % Graph

    pkvec = [];

    for mm=1:length(algs)
        if ~isempty(algs{mm})
            pkvec = union(pkvec,algs{mm}(algs{mm}(:,6) == ii,5));
        end
    end



    for kk=1:1:3 %length(pkvec) % pkvec

        figure(fignum)
        if ii < 3
            if ii == 1
                subplot(2,length(total_pkvec{1}),kk);
            else
                subplot(2,3,kk + 3);
            end
        elseif ii >= 3
            if ii == 3
                subplot(2,3,kk);
            else
                subplot(2,3,kk + 3);
            end
        end

        algplotlist = [];

        for jj=[2 5] % Algorithms (DetRob)

            if ~isempty(algs{jj})
                data = algs{jj}(algs{jj}(:,6) == ii,:);

                pkdata = data(data(:,5) == pkvec(kk),:);
                if ~isempty(pkdata)
                    algplotlist = [algplotlist jj];
                end

                nvec = unique(pkdata(:,4)); % Get list of unique values of n

                avgtime = [];
                maxtime = [];
                mintime = [];
                for ll=1:length(nvec)
                    % Choose all entries with n = nvec(ll)
                    ntime_mat = pkdata(pkdata(:,4) == nvec(ll),:);
                    avgtime(ll) = mean(ntime_mat(:,1));
                    maxtime(ll) = max(ntime_mat(:,1));
                    mintime(ll) = min(ntime_mat(:,1));
                end

                hold on
                errorbar(nvec,avgtime,(avgtime-mintime),(maxtime-avgtime),'o-','color',plotcolors(jj,:))
                hold off
                set(gca, 'YScale', 'log')
                xlabel('Number of nodes')
                if kk == 1
%                     str =
                    ylabel({['\fontsize{16}' graphnames{ii}],'\fontsize{12}Computation time'})
                end
                %                 title(['Graph for ' graphnames(ii) ' graph and pk=' num2str(pkvec(kk))])
                if ii < 3
                    title(['p =' num2str(pkvec(kk))])
                elseif ii >= 3
                    title(['k =' num2str(pkvec(kk))])
                end

                %                 fignum = fignum +1;
            end
        end
        if kk == 1
        legend(algfullnames(algplotlist));
        end

    end

    if mod(ii,2) == 0

        h = get(gcf,'children');
%         haxes = h(2:2:end);
        linkaxes(h,'y');

        fignum = fignum + 1;
    end
end


% LB, rRobustGnl2, rRobupper, rRoblower

algplotlist = [];

for ii=1:4


    pkvec = [];

    for mm=1:length(algs)
        if ~isempty(algs{mm})
            pkvec = union(pkvec,algs{mm}(algs{mm}(:,6) == ii,5));
        end
    end

    for kk=1:1:3 %length(pkvec) % pkvec

        figure(fignum)
        if ii < 3
            if ii == 1
                subplot(2,3,kk);
            else
                subplot(2,3,kk + 3);
            end
        elseif ii >= 3
            if ii == 3
                subplot(2,3,kk);
            else
                subplot(2,3,kk + 3);
            end
        end

        algplotlist = [];

        for jj=[1 4 7 8] % Algorithms

            if ~isempty(algs{jj})
                data = algs{jj}(algs{jj}(:,6) == ii,:);

                pkdata = data(data(:,5) == pkvec(kk),:);
                if ~isempty(pkdata)
                    algplotlist = [algplotlist jj];
                end

                nvec = unique(pkdata(:,4));

                avgtime = [];
                maxtime = [];
                mintime = [];
                for ll=1:length(nvec)
                    ntime_mat = pkdata(pkdata(:,4) == nvec(ll),:);
                    avgtime(ll) = mean(ntime_mat(:,1));
                    maxtime(ll) = max(ntime_mat(:,1));
                    mintime(ll) = min(ntime_mat(:,1));
                end

%                 figure(fignum)
                hold on
                errorbar(nvec,avgtime,(avgtime-mintime),(maxtime-avgtime),'o-','color',plotcolors(jj,:))
                hold off
                set(gca, 'YScale', 'log')
                xlabel('Number of nodes')

                if kk == 1
                ylabel({['\fontsize{16}' graphnames{ii}],'\fontsize{12}Computation time'})
                end
                %                 title(['Graph for ' graphnames(ii) ' graph and pk=' num2str(pkvec(kk))])
                if ii < 3
                    title(['p =' num2str(pkvec(kk))])
                elseif ii >= 3
                    title(['k =' num2str(pkvec(kk))])
                end
            end
        end

        if kk == 1
            legend(algfullnames(algplotlist));
        end
    end

    if mod(ii,2) == 0
        h = get(gcf,'children');
%         haxes = h(2:2:end);
        linkaxes(h,'y');

        fignum = fignum + 1;
    end

end


% Plot the differences between objective values of rRobustGnl2, rRoblower,
% and rRob upper





end
