function plot_robustness(points,num_nodes,fignumber)


if nargin == 3
    figure(fignumber)
else
    figure(1)
end

clf

% Infer additional points

oldpoints = points;

% Plot the values (r,s)-robust ==> (r-1,s+1)-robust
for ii=1:size(oldpoints,1)
    
    r_point = oldpoints(ii,1);
    s_point = oldpoints(ii,2);
    
    for kk = 1:1:min(r_point, num_nodes - s_point)
        if ~ismember([r_point - kk, s_point + kk],points,'rows')
            points = [points; r_point-kk s_point+kk];
        end
    end
end

% Repeat for values (r,s)-robust ==> 1 <= r' <= r and 1 <= s' <= s
oldpoints = points;

for ii=1:size(oldpoints,1) % Iterate through original points
    
    % Plot the values (r,s)-robust ==> 1 <= r' <= r and 1 <= s' <= s
    % robust.
    r_point = oldpoints(ii,1);
    s_point = oldpoints(ii,2);
    
    for rr = r_point:-1:1 % r value
        for ss = s_point:-1:1
            if ~ismember([rr ss],points,'rows')
                points = [points; rr ss];
            end
        end
    end
    
%     % Plot the values (r,s)-robust ==> (r-1,s+1)-robust
%     for kk = 1:1:min(r_point, num_nodes - s_point)
%         if ~ismember([r_point - kk, s_point + kk],points,'rows')
%             points = [points; r_point-kk s_point+kk];
%         end
%     end
    
end

xlimit = ceil(num_nodes/2) + 1;
ylimit = num_nodes + 1;


hold on
for xx=1:ceil(num_nodes/2)
    for yy=1:num_nodes
        if ismember([xx yy],points,'rows')
            o_handle = plot(xx,yy,'o','Color',[0 .6 0],'MarkerSize',10);
        else
            x_handle = plot(xx,yy,'xr','MarkerSize',10);
        end
    end    
end
hold off

xlim([0,xlimit]);
ylim([0,ylimit]);

xlabel('Values of r')
ylabel('Values of s')

legend([o_handle x_handle],'Is robust','Not robust')

end