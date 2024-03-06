clear all
node_num = 1500;  
dim = 2;
coor = rand(dim, node_num);
radius = .06;

% Compute edges
s = [];
t = [];
for i = 1 : node_num
    for j = i : node_num
        dist = sqrt(sum((coor(:,i)-coor(:,j)).^2));
        if dist < radius
            s = [s i];
            t = [t j];
        end
    end
end

G = graph(s, t, 'omitselfloops');
AdjMatrix = adjacency(G);
AdjMatrixFull = full(AdjMatrix);

% Calculate node values for different frequencies and probabilities
tim = 30;
%calf = [1.5 2 4 6 8 10 12 14 16];
calf = [0.1 0.3 .5 1 2 4 6 8 10];
val_frq = size(calf, 2);
MIX = cell(tim, val_frq);

% Assuming you have already defined val_frq, tim, node_num, coor, and calf

% Initialize MIX as a 3D cell array
MIX = cell(tim, val_frq, 11); % 11 for the 11 probability values

for jj = 1 : val_frq
    for ii = 1 : tim
        n = node_num;
        for pr = 0 : .1 : 1
            Z = binornd(1, pr * ones(n, 1));
            Y = 2 * rand(n, 1) * sqrt(3) - sqrt(3);
            frq = calf(jj) * pi;
            X = (1-Z) .* (sin(frq*coor(1,:)') + sin(frq*coor(2,:)')) + Z .* Y;
            
            % Save the result in MIX using three indices
            MIX{ii, jj, round(10 * pr + 1)} = X;
        end
    end
end


%%

% Choose specific indices for time, frequency, and probability
time_idx = 1;
freq_idx = 1;
prob_idx = 8; % This corresponds to 0.5 as values range from 0 to 1 in steps of 0.1

% Extract node values for the specified indices
values_to_plot = MIX{time_idx, freq_idx, prob_idx};

% Plot the graph
figure(1);
p = plot(G, 'XData', coor(1,:), 'YData', coor(2,:));
p.NodeCData = values_to_plot; % Color nodes based on their values
p.MarkerSize = 10;
p.EdgeColor = 'k'; % Set edges to black
colorbar;
title(sprintf('Graph with Node Values for Time %d, Frequency %d, and Probability %.1f', time_idx, freq_idx, (prob_idx-1)*0.1));



%%


% Parameters for PEG
m = 3;
L = 1; % You can adjust this if needed

% Initialize storage for results
Out_PEG_all = zeros(tim, val_frq, 11); % Adjusted to 3D
npdf_all = cell(tim, val_frq, 11); % Adjusted to 3D
vertexPatterns_all = cell(tim, val_frq, 11); % Adjusted to 3D

% Adjacency matrix for the graph
AdjMatrix = adjacency(G);

% Compute PEG for each MIX
for jj = 1 : val_frq
    for h = 1 : tim
        for q = 1 : 11
            % Extract node values for the current indices
            VALS = MIX{h, jj, q};
            
            % Compute PEG for the current node values
            [Out_PEG, npdf, vertexPatterns] = PEG(VALS', AdjMatrix, m, L);
            {h, jj, q}
            % Store results
            Out_PEG_all(h, jj, q) = Out_PEG; % Adjusted indexing
            npdf_all{h, jj, q} = npdf; % Adjusted indexing
            vertexPatterns_all{h, jj, q} = vertexPatterns; % Adjusted indexing
        end
    end
end

% Now, Out_PEG_all, npdf_all, and vertexPatterns_all contain the results for all MIX.


% Now, Out_PEG_all, npdf_all, and vertexPatterns_all contain the results for all MIX.

%%
% Choose specific indices for time, frequency, and probability
time_idx = 11;
freq_idx = 1;
prob_idx = 5; % This corresponds to 0.5 as values range from 0 to 1 in steps of 0.1

% Extract patterns for the chosen MIX across all vertices
patterns_for_MIX = vertexPatterns_all{time_idx, freq_idx, prob_idx};

% Convert patterns to string format for x-tick labels
pattern_strings = cellfun(@(x) strjoin(string(x),'-'), patterns_for_MIX, 'UniformOutput', false);

% Ensure all elements are character vectors
pattern_strings = cellfun(@char, pattern_strings, 'UniformOutput', false);

% Get unique patterns and their counts
[unique_patterns, ~, idx] = unique(pattern_strings);
counts = accumarray(idx, 1);

% Plot the histogram using bar plot
figure(2);
bar(categorical(unique_patterns), counts);
title(['Distribution of patterns for MIX at time ', num2str(time_idx), ', frequency ', num2str(freq_idx), ', and probability ', num2str((prob_idx-1)*0.1)]);
xlabel('Pattern');
ylabel('Frequency');

% Annotate the entropy value
entropy_value = Out_PEG_all(time_idx, freq_idx, prob_idx);
text(0.5, max(counts)*0.9, ['Entropy: ', num2str(entropy_value)], 'FontSize', 10, 'Color', 'red');



%%
%%
% Initialize storage for results
beta_all = zeros(tim, val_frq, 11);
alpha_all = zeros(tim, val_frq, 11);
gamma_all = zeros(tim, val_frq, 11);
delta_all = zeros(tim, val_frq, 11);
tau_all = zeros(tim, val_frq, 11);

% Iterate over all h, jj, and q
for h = 1:tim
    for jj = 1:val_frq
        for q = 1:11
            [beta, alpha, gamma, delta, tau] = transformPatterns(vertexPatterns_all{h, jj, q});
            
            % Store the results
            beta_all(h, jj, q) = beta/node_num;
            alpha_all(h, jj, q) = alpha/node_num;
            gamma_all(h, jj, q) = gamma/node_num;
            delta_all(h, jj, q) = delta/node_num;
            tau_all(h, jj, q) = tau/node_num;
        end
    end
end


%%
% Assuming mean_alpha and std_alpha are both of size 9x11
% Each row corresponds to a different frequency
% Each column corresponds to a different probability

% Calculate mean and standard deviation of alpha values across all time instances
mean_alpha = mean(alpha_all, 1); % Assuming alpha_all is of size tim x val_frq x 11
std_alpha = std(alpha_all, 0, 1);

% Plot results
C = linspecer(val_frq);

figure(1)
for hh = 1:val_frq
    valuem = 0 : .1 : 1;
    errorbar(valuem, squeeze(mean_alpha(1, hh, :)), squeeze(std_alpha(1, hh, :))/2, 'color', C(hh,:), 'LineWidth', 1.5);
    hold on
end

hold off;
xlabel('Probability \it{p}');
ylabel('\alpha');
lgd = legend('0.1','0.3','0.5','1','2','4','6','8','10', 'Location', 'southeast', 'NumColumns', 3);
title(lgd, 'Frequency (Hz)');
xlim([-.01 1.01]);
ylim([min(mean_alpha(:)) - 0.05, max(mean_alpha(:)) + 0.05]); % Adjust ylim based on your data range
set(gca, 'FontSize', 22);
exportgraphics(gca, 'AlphaFrequencyAndP.eps');
exportgraphics(gca, 'figure_2.png');
% %%
% 
% % Assuming mean_entropy and std_entropy are both of size 9x11
% % Each row corresponds to a different frequency
% % Each column corresponds to a different probability
% 
% % Plot results
% C = linspecer(val_frq);
% 
% figure(1)
% for hh = 1:val_frq
%     valuem = 0 : .1 : 1;
%     errorbar(valuem, mean_entropy(hh, :), std_entropy(hh, :)/2, 'color', C(hh,:), 'LineWidth', 1.5);
%     hold on
% end
% 
% hold off;
% xlabel('Probability \it{p}');
% ylabel('Entropy value');
% lgd = legend('1.5','2','4','6','8','10','12','14','16', 'Location', 'southeast', 'NumColumns', 2);
% title(lgd, 'Frequency');
% xlim([-.01 1.01]);
% ylim([min(mean_entropy(:)) - 0.05, max(mean_entropy(:)) + 0.05]); % Adjust ylim based on your data range
% set(gca, 'FontSize', 20);
% exportgraphics(gca, 'PEGfrequencyandp.eps');
% exportgraphics(gca, 'figure_1.png');


%%
% Initialize storage for results
beta_all = zeros(tim, val_frq, 11);
alpha_all = zeros(tim, val_frq, 11);
gamma_all = zeros(tim, val_frq, 11);
delta_all = zeros(tim, val_frq, 11);
tau_all = zeros(tim, val_frq, 11);

% Iterate over all h, jj, and q
for h = 1:tim
    for jj = 1:val_frq
        for q = 1:11
            [beta, alpha, gamma, delta, tau] = transformPatterns(vertexPatterns_all{h, jj, q});
            
            % Store the results
            beta_all(h, jj, q) = beta/node_num;
            alpha_all(h, jj, q) = alpha/node_num;
            gamma_all(h, jj, q) = gamma/node_num;
            delta_all(h, jj, q) = delta/node_num;
            tau_all(h, jj, q) = tau/node_num;
        end
    end
end


%%
% Assuming mean_alpha and std_alpha are both of size 9x11
% Each row corresponds to a different frequency
% Each column corresponds to a different probability

% Calculate mean and standard deviation of alpha values across all time instances
mean_alpha = mean(alpha_all, 1); % Assuming alpha_all is of size tim x val_frq x 11
std_alpha = std(alpha_all, 0, 1);

% Plot results
C = linspecer(val_frq);

figure(1)
for hh = 1:val_frq
    valuem = 0 : .1 : 1;
    errorbar(valuem, squeeze(mean_alpha(1, hh, :)), squeeze(std_alpha(1, hh, :))/2, 'color', C(hh,:), 'LineWidth', 1.5);
    hold on
end

hold off;
xlabel('Probability \it{p}');
ylabel('\alpha');
lgd = legend('0.1','0.3','0.5','1','2','4','6','8','10', 'Location', 'southeast', 'NumColumns', 3);
title(lgd, 'Frequency (Hz)');
xlim([-.01 1.01]);
ylim([min(mean_alpha(:)) - 0.05, max(mean_alpha(:)) + 0.05]); % Adjust ylim based on your data range
set(gca, 'FontSize', 22);
exportgraphics(gca, 'AlphaFrequencyAndP.eps');
exportgraphics(gca, 'figure_2.png');

%%
% Assuming mean_alpha and std_alpha are both of size 9x11
% Each row corresponds to a different frequency
% Each column corresponds to a different probability

% Calculate mean and standard deviation of alpha values across all time instances
mean_beta = mean(beta_all, 1); % Assuming alpha_all is of size tim x val_frq x 11
std_beta = std(beta_all, 0, 1);

% Plot results
C = linspecer(val_frq);

figure(1)
for hh = 1:val_frq
    valuem = 0 : .1 : 1;
    errorbar(valuem, squeeze(mean_beta(1, hh, :)), squeeze(std_beta(1, hh, :))/2, 'color', C(hh,:), 'LineWidth', 1.5);
    hold on
end

hold off;
xlabel('Probability \it{p}');
ylabel('\beta');
lgd = legend('0.1','0.3','0.5','1','2','4','6','8','10', 'Location', 'southeast', 'NumColumns', 3);
title(lgd, 'Frequency (Hz)');
xlim([-.01 1.01]);
ylim([min(mean_beta(:)) - 0.05, max(mean_beta(:)) + 0.05]); % Adjust ylim based on your data range
set(gca, 'FontSize', 22);
exportgraphics(gca, 'AlphaFrequencyAndP.eps');
exportgraphics(gca, 'figure_3.png');


%%

%%
% Assuming mean_Out_PEG and std_Out_PEG are both of size 9x11
% Each row corresponds to a different frequency
% Each column corresponds to a different probability

% Calculate mean and standard deviation of Out_PEG values across all time instances
mean_Out_PEG = mean(Out_PEG_all, 1); % Assuming Out_PEG_all is of size tim x val_frq x 11
std_Out_PEG = std(Out_PEG_all, 0, 1);

% Plot results
C = linspecer(val_frq);

figure(2) % Changed the figure number to 2 to avoid overwriting previous figure
for hh = 1:val_frq
    valuem = 0 : .1 : 1;
    errorbar(valuem, squeeze(mean_Out_PEG(1, hh, :)), squeeze(std_Out_PEG(1, hh, :))/2, 'color', C(hh,:), 'LineWidth', 1.5);
    hold on
end

hold off;
xlabel('Probability \it{p}');
ylabel('Entropy value'); % Changed y-label to 'Out_PEG'
lgd = legend('0.1','0.3','0.5','1','2','4','6','8','10', 'Location', 'northeast', 'NumColumns', 3);
title(lgd, 'Frequency (Hz)');
xlim([-.01 1.01]);
ylim([min(mean_Out_PEG(:)) - 0.05, max(mean_Out_PEG(:)) + 0.05]); % Adjust ylim based on your data range
set(gca, 'FontSize', 22);
exportgraphics(gca, 'Out_PEG_FrequencyAndP.eps'); % Updated the filename
exportgraphics(gca, 'figure_4.png'); % Updated the figure number


