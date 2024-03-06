function [beta, alpha, gamma, delta, tau] = transformPatterns(patterns)
    % Initialize counts
    count_beta = 0; % Count for patterns 1,2,3 and 3,2,1
    count_alpha = 0; % Count for patterns 1,3,2 + 2,3,1 + 2,1,3 + 3,1,2
    count_gamma = 0; % Count for patterns 2,3,1 + 2,1,3 - 1,3,2 - 3,1,2
    count_delta = 0; % Count for patterns 1,3,2 + 1,2,3 - 2,3,1 - 3,1,2
    
    % Loop through each pattern and update counts
    for i = 1:length(patterns)
        pattern = patterns{i};
        
       if isequal(pattern, [1,2,3])
            count_beta = count_beta + 1;
        elseif isequal(pattern, [3,2,1])
            count_beta = count_beta - 1;
        elseif isequal(pattern, [1,3,2]) || isequal(pattern, [2,3,1]) || isequal(pattern, [2,1,3]) || isequal(pattern, [3,1,2])
            count_alpha = count_alpha + 1;
        end
        
        if isequal(pattern, [2,3,1]) || isequal(pattern, [2,1,3])
            count_gamma = count_gamma + 1;
        elseif isequal(pattern, [1,3,2]) || isequal(pattern, [3,1,2])
            count_gamma = count_gamma - 1;
        end
        
        if isequal(pattern, [1,3,2]) || isequal(pattern, [1,2,3])
            count_delta = count_delta + 1;
        elseif isequal(pattern, [2,3,1]) || isequal(pattern, [3,1,2])
            count_delta = count_delta - 1;
        end
    end
    
    % Compute the final values
    beta = count_beta;
    alpha = count_alpha;
    gamma = count_gamma;
    delta = count_delta;
    tau = 1/3 - alpha;
end
