function outstruct = rRobB2rpt(args)

% Use the same args as rRobBench2, but with the additional field
% args.repeat : integer of how many times to repeat rRobBenchmark2

% Make a directory to store the results

dirname = ['./data/' date]; % Save dirname to a char because sometimes results take multiple days to compute

if exist(dirname) ~= 7
    mkdir(dirname);
end

save([dirname '/args'],'args');

for ii=1:1:args.repeat
    
    results = struct();
    
    results = rRobBenchmark2(args);
    
    save([dirname '/results' num2str(ii)],'results');
    
end


end