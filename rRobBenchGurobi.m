function outstruct = rRobBenchGurobi(args)

% args.testedMatrices : Matrix of testedMatrices info
% args.testWarm : Set to 0 to skip the warm start method, otherwise leave
% blank or set to 1

testedMatrices = args.testedMatrices;

if (isfield(args,'testWarm') && args.testWarm == 0)
    testWarm = 0;
else
    testWarm = 1;
end

% Close parpool, test Gurobi implementations on all saved matrices

pool = gcp('nocreate')

if ~isempty(pool)
    delete(pool)
end

% Gurobi implementation

time1 = 0;
time2 = 0;

finalgurobi.Erdos = [];
finalgurobi.complete = [];
finalgurobi.kundir = [];
finalgurobi.kdir = [];
finalgurobi.randdir = [];
finalgurobi.kinrand = [];
finalgurobi.koutrand = [];

if testWarm == 1
    finalgurobiwm.Erdos = [];
    finalgurobiwm.complete = [];
    finalgurobiwm.kundir = [];
    finalgurobiwm.kdir = [];
    finalgurobiwm.randdir = [];
    finalgurobiwm.kinrand = [];
    finalgurobiwm.koutrand = [];
end

totmat = size(testedMatrices,1);
disp([newline 'Total matrices to test:' ' ' num2str(totmat) newline])

for ii=1:1:size(testedMatrices,1)
    
    
    tempcell = testedMatrices(ii,:);
    
    L = tempcell{1};
    n = tempcell{3};
    
    if strcmp(tempcell{2},'Erdos')
        p = tempcell{4};
        
        tic
        gstruct = rRobustGnl_gurobi(struct('L',L));
        time1 = toc;
        
        finalgurobi.Erdos = [finalgurobi.Erdos; [time1 gstruct.maxr n p]];
        
        if testWarm == 1
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobiwm.Erdos = [finalgurobiwm.Erdos; [time2 gstructwm.maxr n p]];
        end
        
        
    elseif strcmp(tempcell{2},'complete')
        
        tic
        gstruct = rRobustGnl_gurobi(struct('L',L));
        time1 = toc;
        
        finalgurobi.complete = [finalgurobi.complete; [time1 gstruct.maxr n]];
        
        if testWarm == 1
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobiwm.complete = [finalgurobiwm.complete; [time2 gstructwm.maxr n]];
        end
        
        
    elseif strcmp(tempcell{2},'kundir')
        k = tempcell{4};
        
        tic
        gstruct = rRobustGnl_gurobi(struct('L',L));
        time1 = toc;
        
        finalgurobi.kundir = [finalgurobi.kundir; [time1 gstruct.maxr n k]];
        
        if testWarm == 1
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobiwm.kundir = [finalgurobiwm.kundir; [time2 gstructwm.maxr n k]];
        end
        
    elseif strcmp(tempcell{2},'kdir')
        k = tempcell{4};
        
        tic
        gstruct = rRobustGnl_gurobi(struct('L',L));
        time1 = toc;
        
        finalgurobi.kdir = [finalgurobi.kdir; [time1 gstruct.maxr n k]];
        
        if testWarm == 1
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobiwm.kdir = [finalgurobiwm.kdir; [time2 gstructwm.maxr n k]];
        end
        
    elseif strcmp(tempcell{2},'randdir')
        p = tempcell{4};
        
        tic
        gstruct = rRobustGnl_gurobi(struct('L',L));
        time1 = toc;
        
        finalgurobi.randdir = [finalgurobi.randdir; [time1 gstruct.maxr n p]];
        
        if testWarm == 1
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobiwm.randdir = [finalgurobiwm.randdir; [time2 gstructwm.maxr n p]];
        end
        
    elseif strcmp(tempcell{2},'kinrand')
        k = tempcell{4};
        
        tic
        gstruct = rRobustGnl_gurobi(struct('L',L));
        time1 = toc;
        
        finalgurobi.kinrand = [finalgurobi.kinrand; [time1 gstruct.maxr n k]];
        
        if testWarm == 1
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobiwm.kinrand = [finalgurobiwm.kinrand; [time2 gstructwm.maxr n k]];
        end 
        
    elseif strcmp(tempcell{2},'koutrand')
        k = tempcell{4};
        
        tic
        gstruct = rRobustGnl_gurobi(struct('L',L));
        time1 = toc;
        
        finalgurobi.koutrand = [finalgurobi.koutrand; [time1 gstruct.maxr n k]];
        
        if testWarm == 1
            tic
            gstructwm = rRobustWm_gurobi(struct('L',L));
            time2 = toc;
            
            finalgurobiwm.koutrand = [finalgurobiwm.koutrand; [time2 gstructwm.maxr n k]];
        end
    end
    
    
    if mod(ii,10) == 0
        disp([num2str(ii) ' / ' num2str(totmat) ' matrices analyzed'])
    end
    
end

outstruct.finalgurobi = finalgurobi;

if testWarm == 1
    outstruct.finalgurobiwm = finalgurobiwm;
end

disp('Gurobi analysis done')



end