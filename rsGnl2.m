function outstruct = rsGnl2(args,varargin)

% Determines the max values of r and s for which a digraph is (r,s)-robust using intlinprog. Also measures
%   the computation time.
% This formulation first solves for rho(D) using the "infinity norm" method.
%   It then solves for s using the ~(AvBvC) method with r fixed at rho(D).
%
% args must be a struct with the following fields:
%   args.L : Graph Laplacian
%   args.print : (Optional) Controls whether output is displayed or not. Set
%       to 1 to print output of intlinprog and final results to the screen. Set
%       to 0 for no printed output. Default value is 0.
%   args.init : initial feasible point for the optimization problem
%   args.MaxTime : Maximum time limit for algorithm to run. Do not define
%                   this field for default value of 1e20

% varargin : set as a value of r to only check s for that particular value
%            of r

L = args.L;

if size(L,1) ~= size(L,2)
    error('Matrix is not square')
end

time = 0;
time2 = 0;

first_timer = tic;

n = size(L,1);
nr2 = ceil(n/2);
exitflag2 = 0;
exitflagD = 0;
indegflag = 0; % Equal to 1 if the min. in-degree is high enough to prove that s == n.
not_solved_before_MaxTime_bool = 0;

% Approx bounds for r optimization
approx_lb_r = -1;
approx_ub_r = Inf;

% Approx bounds for s optimization
approx_lb_s = -1;
approx_ub_s = Inf;



    function stop = TimeLimit(x,optimValues,state)
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

% This is a separate function in case I want to save separate bounds
% for r and s optimization cases
    function stop = TimeLimit2(x,optimValues,state)
        if toc(first_timer) > args.MaxTime
            stop = true;
            if ~isequal(state,'done')
                not_solved_before_MaxTime_bool = 1;
                if strcmp(optimValues.phase,'rootlp')
                    % See MATLAB documentation; the lb and ub aren't
                    % necessarily set in rootlp phase
                    approx_lb_s = -1;
                    approx_ub_s = Inf;
                else
                    approx_lb_s = optimValues.lowerbound;
                    approx_ub_s = optimValues.fval;
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

if size(varargin) == 0
    
    % Solve for rho(D)
    
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
    
    
    
    % Display or Suppress iterations, based on args.print
    if ~isfield(args,'print') || args.print == 0
        options.Display = 'off';
    end
    
    
    
    [x,fval,exitflag,output] = intlinprog(c,2:2*n+1,A,b,[],[],lb,ub,options);
    % [x,fval,exitflag,output] = linprog(c,[],[],Aeq,beq,lb,ub,options);
    % [x,fval,exitflag,output] = linprog(c,A,b,[],[],lb,ub,options);
    
    
    r = x(1);
    
elseif size(varargin) == 1
    r = varargin{1};
else
    error('Too many input arguments')
end

time = toc(first_timer);

if nargin == 1
    if exitflag == 2 || 0
        not_solved_before_MaxTime_bool = 1;
    end
end

second_timer = tic;

% Solve for max s for which graph is (rho(D),s)-robust

% Only run this portion if MaxTime has
% not been exceeded and the graph is more than 0-robust
if not_solved_before_MaxTime_bool == 0 && r > 0
    
    degmin = min(diag(L)); % Get minimum in-degree of the graph
    
    
    % See Property 5.23 in LeBlanc 2013, PhD Thesis. This tests
    if degmin >= floor(n/2) + r - 1
        disp('Minimum in-degree greater than (floor(n/2) + r - 1); therefore s = n. No optimization necessary.')
        outstruct.s = n;
        outstruct.s_ub = n;
        indegflag = 1;
    else
        
        % Constraint matrices
        AC = [0 -ones(1,n) zeros(1,3*n);
            0 ones(1,n) zeros(1,3*n);
            0 zeros(1,n) -ones(1,n) zeros(1,2*n);
            0 zeros(1,n) ones(1,n) zeros(1,2*n);
            zeros(n,1) eye(n) eye(n) zeros(n,n) zeros(n,n)];
        
        % Unnecessary constraints
        % BC = [-ones(n,1) zeros(n,1) L zeros(n,3*n);
        %     -ones(n,1) zeros(n,1) zeros(n) L zeros(n,2*n)];
        
        CC = [0 -ones(1,n) zeros(1,n) ones(1,n) zeros(1,n);
            0 zeros(1,n) -ones(1,n) zeros(1,n) ones(1,n);
            -1 zeros(1,n) zeros(1,n) ones(1,n) ones(1,n)];
        
        DC = [zeros(n,1) L zeros(n) -(n-1)*eye(n) zeros(n);
            zeros(n,1) zeros(n) L  zeros(n) -(n-1)*eye(n)];
        
        EX = [0 0 zeros(1,2*n) -ones(1,2*n)]; % testing only
        
        
        A2 = [AC; CC; DC];
        
        bAC = [-1; (n-1); -1; (n-1); ones(n,1)];
        %         bBC = [zeros(2*n,1)];
        bCC = [-1;-1;0];
        bDC = [(r-1)*ones(2*n,1)];
        
        b2 = [bAC; bCC; bDC];
        
        lb2 = [0; zeros(4*n,1)];
        ub2 = [Inf; ones(4*n,1)];
        
        c2 = [1; zeros(4*n,1)];
        
        % Dual formulation
        
        lbd = zeros(size(A2,1));
        
        if isfield(args,'MaxTime')
            options.MaxTime = args.MaxTime;
            options.OutputFcn = @TimeLimit2;
        else
            options.MaxTime = 1e20;
        end
        
        % First test to see if the dual of the convex relaxation is
        % unbounded above. If it is, then the integer formulation is
        % infeasible, implying s = n.
        
        % Dual problem:
        % THIS IS WRONG. Do not use. kept here simply for the record.
        %         [xD, fvalD, exitflagD, outputD] = linprog(b2,[],[],A2',-c2,lbd,[],options);
        %
        %         if exitflagD == -3 || (exitflagD == 1 && fvalD >= n)
        %             time2 = toc;
        %             disp('Problem infeasible (dual opt val >= n); s == n')
        %             outstruct.s = n;
        %         else
        [x2,fval2,exitflag2,output2] = intlinprog(c2,2:4*n+1,A2,b2,[],[],lb2,ub2,options);
        % [x,fval,exitflag,output] = linprog(c,[],[],Aeq,beq,lb,ub,options);
        % [x,fval,exitflag,output] = linprog(c,A,b,[],[],lb,ub,options);
        
        
        % If the problem is infeasible, then s == n.
        if exitflag2 == -2
            disp('Problem infeasible (primal of s problem); s == n')
            outstruct.s = n;
            outstruct.s_ub = n;
        end
        %         end
    end
else % This option reached if r <= 0 or MaxTime violated.
    if r <= 0
        disp('Graph is 0-robust.')
        outstruct.s = 0;
        outstruct.s_ub = 0;
    end
end

time2 = toc(second_timer);

if exitflag2 == 2 || 0
    not_solved_before_MaxTime_bool = 1;
end


% Only put this data in outstruct if MaxTime hasn't been violated

if not_solved_before_MaxTime_bool == 0
    
    outstruct.success_before_MaxTime_bool = 1;
    
    outstruct.r = r;
    outstruct.r_ub = r;
    
    if size(varargin) == 0
        outstruct.b1 = x(2:n+1);
        outstruct.b2 = x(n+2:2*n+1);
        
        outstruct.exitflag = exitflag;
        outstruct.output = output;
        outstruct.time = time;
        outstruct.totaltime = time + time2;
        
        if exitflag == 2 || exitflag == 0
            outstruct.approx_lb = approx_lb_r;
            outstruct.approx_ub = approx_ub_r;
        end
    end
    
    if exitflag2 == 1
        outstruct.s = x2(1);
        outstruct.s_ub = outstruct.s;
        
        outstruct.two.b1 = x2(2:n+1);
        outstruct.two.b2 = x2(n+2:2*n+1);
        outstruct.two.y1 = x2(2*n+2:3*n+1);
        outstruct.two.y2 = x2(3*n+2:end);
        
        outstruct.two.exitflag = exitflag2;
        outstruct.two.output = output2;
        outstruct.two.time = time2;
        % elseif exitflag2 == 2 || exitflag2 == 0
        
    elseif exitflag2 == -2
        % Problem infeasible
        disp('Again, problem infeasible (primal of s problem); s == n')
    elseif indegflag == 0
        % !!!! THIS SECTION IS QUESTIONABLE. If there are problems, check
        %   here first.
        disp(['No data for second optimization program because exitflag2 = ' num2str(exitflag2)])
        if ~isfield(outstruct,'s')
            % Since this will only be reached if r > 0, then 
            outstruct.s = 1;
            outstruct.s_ub = approx_ub_s;
        end
    end
    
    
    if exitflagD ~= 0 && exitflag ~= -3
        outstruct.dual.exitflag = exitflagD;
        outstruct.dual.output = outputD;
        outstruct.dual.fval = fvalD;
        outstruct.dual.x = xD;
    end
    
    if isfield(args,'print') && args.print == 1
        disp(['Max r:' ' ' num2str(r)])
        disp(['Max s:' ' ' num2str(outstruct.s)])
        disp(['Total opt time: ' ' ' num2str(time + time2)])
        disp(['First opt time: ' ' ' num2str(time)])
        if r > 0
            disp(['Second opt time: ' ' ' num2str(time2)])
        end
        
        if exitflag2 == 1
            disp('One possible minimal S1 / S2 pair:')
            disp(['S1:' ' ' num2str(round(outstruct.two.b1'))])
            disp(['S2: ' ' ' num2str(round(outstruct.two.b2'))])
        else
            disp('One possible minimal S1 / S2 pair:')
            disp(['S1:' ' ' num2str(round(outstruct.b1'))])
            disp(['S2:' ' ' num2str(round(outstruct.b2'))])
        end
    end
    
else % If MaxTime was violated
    
    outstruct.success_before_MaxTime_bool = 0;
    
    outstruct.r = approx_lb_r;
    outstruct.r_ub = approx_ub_r;
    outstruct.s_ub = approx_ub_s;
    outstruct.s = approx_lb_s;
    outstruct.totaltime = args.MaxTime;
    
end

end