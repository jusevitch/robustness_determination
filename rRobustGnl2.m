function outstruct = rRobustGnl2(args)

% General integer linear program r-robustness determination function.
% Determine the r-robustness of a digraph using intlinprog. Also measures
%   the computation time.
% *** This formulation uses less optimization constraints than the prior function.
% args must be a struct with the following fields:
%   args.L : Graph Laplacian
%   args.print : (Optional) Controls whether output is displayed or not. Set
%       to 1 to print output of intlinprog and final results to the screen. Set
%       to 0 for no printed output. Default value is 0.
%   args.init : initial feasible point for the optimization problem
%   args.MaxTime : Maximum time limit for algorithm to run. Do not define
%                   this field for default value of 1e20

L = args.L;

if size(L,1) ~= size(L,2)
    error('Matrix is not square')
end

time = 0;

first_timer = tic;

n = size(L,1);
nr2 = ceil(n/2);

not_solved_before_MaxTime_bool = 0;

% Approx bounds for r optimization
approx_lb_r = -1;
approx_ub_r = Inf;

e1 = [1; zeros(n-1,1)];

A = [-ones(n,1) L zeros(n);
    -ones(n,1) zeros(n) L;
    0 -ones(1,n) zeros(1,n);
    0 ones(1,n) zeros(1,n);
    0 zeros(1,n) -ones(1,n);
    0 zeros(1,n) ones(1,n);
    zeros(n,1) eye(n) eye(n)];

b = [zeros(2*n,1); -1; (n-1); -1; (n-1); ones(n,1)];

lb = [0; zeros(2*n,1)];
ub = [ceil(n/2); ones(2*n,1)]; % Seems to run faster with the upper bound


c = [1; zeros(2*n,1)];

    function stop = TimeLimit(x,optimValues,state)
        %         toc(first_timer)
        if toc(first_timer) > args.MaxTime
            stop = true;
            if ~isequal(state,'done')
                not_solved_before_MaxTime_bool = 1;
                if strcmp(optimValues.phase,'rootlp')
                    % See MATLAB documentation; the lb and ub aren't
                    % necessarily set in rootlp phase
                    approx_lb_r = -1;
                    approx_ub_r = Inf;
                else
                    approx_lb_r = optimValues.lowerbound;
                    approx_ub_r = optimValues.fval;
                end
            end
        else
            stop = false;
        end
    end

% Set Integer tolerance to the lowest possible value
% Change the Heuristics to basic or intermediate to possibly save time
options = optimoptions('intlinprog','IntegerTolerance',1e-6,'MaxNodes',2^(3*n),'Heuristics','advanced');
if isfield(args,'MaxTime')
    options.MaxTime = args.MaxTime;
    options.OutputFcn = @TimeLimit;
else
    options.MaxTime = 1e20;
end




% Display or Suppress iterations, based on args.print
if ~isfield(args,'print') || args.print == 0
    options.Display = 'off';
end



[x,fval,exitflag,output] = intlinprog(c,2:2*n+1,A,b,[],[],lb,ub,options);
% [x,fval,exitflag,output] = linprog(c,[],[],Aeq,beq,lb,ub,options);
% [x,fval,exitflag,output] = linprog(c,A,b,[],[],lb,ub,options);

time = toc(first_timer);

if exitflag == 2 || 0
    not_solved_before_MaxTime_bool = 1;
end

if not_solved_before_MaxTime_bool == 0 % If solution was found before MaxTime was reached
    
    outstruct.success_before_MaxTime_bool = 1;
    
    outstruct.r = x(1);
    outstruct.r_ub = outstruct.r;
    outstruct.b1 = x(2:n+1);
    outstruct.b2 = x(n+2:2*n+1);
    % outstruct.b0 = x(2*n+2:end);
    outstruct.exitflag = exitflag;
    outstruct.output = output;
    outstruct.time = time;
    
    if isfield(args,'print') && args.print == 1
        disp('Max r: ')
        x(1)
        disp('Computation time:')
        time
        
        disp('One possible minimal S1 / S2 pair:')
        disp('S1:')
        outstruct.b1'
        disp('S2: ')
        outstruct.b2'
        %     disp('S0: ')
        %     outstruct.b0'
    end
    
else % If solution was not found before MaxTime was reached
    
    outstruct.success_before_MaxTime_bool = 0;
    
    outstruct.r = approx_lb_r;
    outstruct.r_ub = approx_ub_r;
    outstruct.s = -1;
    outstruct.time = args.MaxTime;
    
end

end