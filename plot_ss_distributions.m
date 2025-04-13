output_dir = fullfile(pwd, 'figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% === Plot: Het TFP Model ===
load('steady_state_firmdata_TFP_final.mat', 'z_vec', 'k_firms', 'n_firms', 'y_firms');

% Compute output concentration: share of output from top 10% most productive firms
N = length(y_firms);
[~, sorted_idx] = sort(y_firms, 'descend');
top10_idx = sorted_idx(1:ceil(0.10 * N));
output_share_top10 = sum(y_firms(top10_idx)) / sum(y_firms);
fprintf('Top 10%% of firms account for %.1f%% of total output in the het-TFP model.\n', 100 * output_share_top10);

% === Capital concentration (top 10%)
[~, sorted_k_idx] = sort(k_firms, 'descend');
top10_k_idx = sorted_k_idx(1:ceil(0.10 * N));
capital_share_top10 = sum(k_firms(top10_k_idx)) / sum(k_firms);
fprintf('Top 10%% of firms hold %.1f%% of total capital in the het-TFP model.\n', 100 * capital_share_top10);

% === Labor concentration (top 10%)
[~, sorted_n_idx] = sort(n_firms, 'descend');
top10_n_idx = sorted_n_idx(1:ceil(0.10 * N));
labor_share_top10 = sum(n_firms(top10_n_idx)) / sum(n_firms);
fprintf('Top 10%% of firms employ %.1f%% of total labor in the het-TFP model.\n', 100 * labor_share_top10);

figure('Name', 'Steady State Distributions (Het TFP)', ...
       'Color', 'w', 'Position', [100, 100, 1000, 800]);

subplot(2,2,1);  plot_with_density(z_vec, 15, 'Productivity Distribution', 'z_i', [0.2 0.6 1]);
subplot(2,2,2);  plot_with_density(k_firms, 15, 'Capital Distribution', 'k_i', [0.3 0.7 0.4]);
subplot(2,2,3);  plot_with_density(n_firms, 15, 'Labor Distribution', 'n_i', [0.7 0.4 0.4]);
subplot(2,2,4);  plot_with_density(y_firms, 15, 'Output Distribution', 'y_i', [0.6 0.5 0.9]);

% === Export to PDF
if ~exist('figures', 'dir'); mkdir('figures'); end
exportgraphics(gcf, fullfile(output_dir, 'het_TFP_distributions.pdf'), 'ContentType', 'vector');

%% === Plot: Het Xi Model ===
load('steady_state_firmdata_xi_final.mat', 'xi_vec', 'k_firms', 'n_firms', 'y_firms');

% === Compute output concentration (top 10%)
N = length(y_firms);
[~, sorted_y_idx] = sort(y_firms, 'descend');
top10_y_idx = sorted_y_idx(1:ceil(0.10 * N));
output_share_top10 = sum(y_firms(top10_y_idx)) / sum(y_firms);
fprintf('Top 10%% of firms account for %.1f%% of total output in the het-\\xi model.\n', 100 * output_share_top10);

% === Capital concentration (top 10%)
[~, sorted_k_idx] = sort(k_firms, 'descend');
top10_k_idx = sorted_k_idx(1:ceil(0.10 * N));
capital_share_top10 = sum(k_firms(top10_k_idx)) / sum(k_firms);
fprintf('Top 10%% of firms hold %.1f%% of total capital in the het-\\xi model.\n', 100 * capital_share_top10);

% === Labor concentration (top 10%)
[~, sorted_n_idx] = sort(n_firms, 'descend');
top10_n_idx = sorted_n_idx(1:ceil(0.10 * N));
labor_share_top10 = sum(n_firms(top10_n_idx)) / sum(n_firms);
fprintf('Top 10%% of firms employ %.1f%% of total labor in the het-\\xi model.\n', 100 * labor_share_top10);

figure('Name', 'Steady State Distributions (Het Xi)', ...
       'Color', 'w', 'Position', [150, 150, 1000, 800]);

subplot(2,2,1);  plot_with_density(xi_vec, 15, 'Financial Conditions Distribution', '\xi_i', [0.5 0.2 0.6]);
subplot(2,2,2);  plot_with_density(k_firms, 15, 'Capital Distribution', 'k_i', [0.3 0.7 0.4]);
subplot(2,2,3);  plot_with_density(n_firms, 15, 'Labor Distribution', 'n_i', [0.7 0.4 0.4]);
subplot(2,2,4);  plot_with_density(y_firms, 15, 'Output Distribution', 'y_i', [0.6 0.5 0.9]);

% === Export to PDF
if ~exist('figures', 'dir'); mkdir('figures'); end
exportgraphics(gcf, 'figures/het_xi_distributions.pdf', 'ContentType', 'vector');

%% === Define Helper Plot Function
function plot_with_density(data, bins, titleStr, xlab, color)
    % Plot histogram normalized to show firm shares
    h = histogram(data, bins, 'Normalization', 'probability', ...
                  'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.4); 
    hold on;

    % Kernel density estimate
    [f, xi] = ksdensity(data);
    binWidth = mean(diff(h.BinEdges));
    f_scaled = f * binWidth;  % Rescale density to match histogram

    plot(xi, f_scaled, 'k-', 'LineWidth', 1.5);

    % Compute and plot mean
    m = mean(data);
    xline(m, 'r--', 'LineWidth', 1.5);
    text(m, max(f_scaled)*0.9, sprintf('Mean: %.2f', m), ...
         'Color', 'r', 'FontSize', 8, 'HorizontalAlignment', 'left', 'Rotation', 90);

    % Lognormal fit
    if strcmp(xlab, 'z_i') || strcmp(xlab, '\xi_i')
        pd = fitdist(data, 'Lognormal');
        y_fit = pdf(pd, xi) * binWidth;
        plot(xi, y_fit, '--', 'Color', [0 0.5 0], 'LineWidth', 1.2);
        legend('Histogram', 'Kernel Density', 'Lognormal Fit', 'Mean', 'Location', 'best');
    else
        legend('Histogram', 'Kernel Density', 'Mean', 'Location', 'best');
    end

    fprintf('Mean of %s: %.4f\n', xlab, m);

    title(titleStr); 
    xlabel(xlab); 
    ylabel('Fraction of Firms'); 
    grid on;
end


