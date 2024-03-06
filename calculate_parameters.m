function [alpha, beta, gamma, delta, tau] = calculate_parameters(vec)
    % Check if the input is a 1x6 vector
    if size(vec, 1) ~= 1 || size(vec, 2) ~= 6
        error('Input must be a 1x6 vector');
    end
    
    % Calculate alpha
    alpha = vec(5) + vec(3) + vec(4) + vec(2);
    
    % Calculate beta
    beta = vec(6) - vec(1);
    
    % Calculate gamma
    gamma = vec(3) + vec(4) - vec(5) - vec(2);
    
    % Calculate delta
    delta = vec(5) + vec(4) - vec(3) - vec(2);
    
    % Calculate tau
    tau = (2/3) - alpha;
end
