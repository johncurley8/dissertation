% Load IRFs
clear oo_ irfs_het
load('./Baseline RBC/IRFs_JQ.mat');       irfs_JQ     = oo_.irfs;
load('IRFs_het_xi_final.mat');                  irfs_het_xi = oo_.irfs;
load('IRFs_het_TFP.mat');                irfs_het_TFP = oo_.irfs;

% Create 'figures' folder if it doesn't exist
fig_dir = fullfile(pwd, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Aesthetic Settings
dark_green = [0 0.5 0];
dark_blue = [0 0 0.7];
style_jq   = {'r', 'LineWidth', 1.4};
style_het_xi = {'--', 'Color', dark_green, 'LineWidth', 1.4};
style_het_TFP = {'-.', 'Color', dark_blue, 'LineWidth', 1.4};
timeRange  = 0:30;

y_limits = [-1 1; -1.2 1.2; -1 1; -4 4; -2.5 2.5; -2.0 2.0];

% === Figure: Combined Response to Financial Shock
fig_comb = figure('Name','Response to Financial Shock (Combined)', ...
                  'Color','w','Position',[100 100 1000 600]);

irf_vars = {'yhat_eps_xi', 'nhat_eps_xi', 'chat_eps_xi', ...
            'ihat_eps_xi', 'byhat_eps_xi', 'dyhat_eps_xi'};

var_titles = {'Output', 'Hours', 'Consumption', ...
              'Investment', 'Debt rep./Y', 'Equity Payout/Y'};

label_font = 10;
title_font = 11;

for i = 1:6
    subplot(2,3,i)
    plot(-irfs_JQ.(irf_vars{i})*100, style_jq{:}); hold on
    plot(-irfs_het_xi.(irf_vars{i})*100, style_het_xi{:})
    plot(-irfs_het_TFP.(irf_vars{i})*100, style_het_TFP{:})
    
    axis([0 30 y_limits(i,:)])
    title(var_titles{i}, 'FontSize', title_font)
    xlabel('Quarters', 'FontSize', label_font)
    ylabel('Percent', 'FontSize', label_font)
    grid on; box off

    if i == 1
        legend('JQ (rep. firm; constraint unch.)', 'Heterogeneous firms (\xi)', ...
               'Heterogeneous firms (z)', ...
               'Location','best', 'FontSize', 10);
    end
end

% Export as vector PDF
exportgraphics(fig_comb, 'figures/IRF_finshock_combined.pdf', 'ContentType', 'vector');

% === Figure: Combined Response to TFP Shock
fig_comb_tfp = figure('Name','Response to TFP Shock (Combined)', ...
                      'Color','w','Position',[100 100 1000 600]);

irf_vars = {'yhat_eps_z', 'nhat_eps_z', 'chat_eps_z', ...
            'ihat_eps_z', 'byhat_eps_z', 'dyhat_eps_z'};

for i = 1:6
    subplot(2,3,i)
    plot(-irfs_JQ.(irf_vars{i})*100, style_jq{:}); hold on
    plot(-irfs_het_xi.(irf_vars{i})*100, style_het_xi{:})
    plot(-irfs_het_TFP.(irf_vars{i})*100, style_het_TFP{:})
    
    axis([0 30 y_limits(i,:)])
    title(var_titles{i}, 'FontSize', title_font)
    xlabel('Quarters', 'FontSize', label_font)
    ylabel('Percent', 'FontSize', label_font)
    grid on; box off

    if i == 1
        legend('JQ (rep. firm; constraint unch.)', 'Heterogeneous firms (\xi)', ...
               'Heterogeneous firms (z)', ...
               'Location','best', 'FontSize', 10);
    end
end

% Export as vector PDF
exportgraphics(fig_comb_tfp, 'figures/IRF_tfpshock_combined.pdf', 'ContentType', 'vector');

%}


% === Setup
N = 250;               % Number of firms
T = 30;               % IRF horizon
time = 0:T-1;

var_base = {'yhat', 'nhat', 'ihat', 'byhat', 'dyhat'};
titles   = {'Output', 'Labor', 'Investment', 'Debt rep/yi', 'Equity payout/yi'};

% === Style settings
style_tfp = {'--', 'Color', dark_green, 'LineWidth', 1.4};
style_fin = {'-', 'Color', dark_blue, 'LineWidth', 1.4};

% === Create figure
fig_mean_compare = figure('Name','Mean Firm IRFs to Idiosyncratic Shocks (TFP vs Financial)', ...
                          'Color','w', 'Position', [100 100 1000 600]);

for v = 1:length(var_base)
    mean_tfp = zeros(1, T);
    mean_fin = zeros(1, T);
    varname_base = var_base{v};

    for i = 1:N
        tfp_irf = sprintf('%s%d_eps_z%d', varname_base, i, i);
        fin_irf = sprintf('%s%d_eps_xi%d', varname_base, i, i);

        if isfield(irfs_het_TFP, tfp_irf)
            mean_tfp = mean_tfp + irfs_het_TFP.(tfp_irf)(1:T);
        else
            warning('Missing TFP field: %s', tfp_irf);
        end

        if isfield(irfs_het_xi, fin_irf)
            mean_fin = mean_fin + irfs_het_xi.(fin_irf)(1:T);
        else
            warning('Missing Financial field: %s', fin_irf);
        end
    end

    mean_tfp = mean_tfp / N;
    mean_fin = mean_fin / N;

    subplot(2, 3, v);
    plot(time, -mean_tfp * 100, style_tfp{:}); hold on;
    plot(time, -mean_fin * 100, style_fin{:});
    title(['Mean Firm ', titles{v}, ' Response']);
    xlabel('Quarters'); ylabel('Percent');
    grid on; box off;
    
    y_max = 1.1 * max(abs([mean_tfp, mean_fin] * 100), [], 'all');
    axis([0 T -y_max y_max]);

    if v == 1
        legend('TFP Shock', 'Financial Shock', 'Location', 'best', 'FontSize', 9);
    end
end
exportgraphics(gcf, 'figures/IRF_mean_idio_compare.pdf', 'ContentType', 'vector');
%}


% === TFP Shock: Mean + Percentile IRFs ===

var_base = {'yhat', 'nhat', 'ihat', 'byhat', 'dyhat'};
titles   = {'Output', 'Labor', 'Investment', 'Debt rep/yi', 'Equity payout/yi'};

fig_tfp = figure('Name','TFP Shock: Mean + Percentiles', ...
                 'Color','w', 'Position', [100 100 1000 600]);

for v = 1:length(var_base)
    tfp_mat = nan(N, T);
    varname_base = var_base{v};

    for i = 1:N
        tfp_irf = sprintf('%s%d_eps_z%d', varname_base, i, i);
        if isfield(irfs_het_TFP, tfp_irf)
            tfp_mat(i, :) = irfs_het_TFP.(tfp_irf)(1:T);
        end
    end

    mean_tfp = mean(tfp_mat, 1, 'omitnan');
    pct25 = prctile(tfp_mat, 25, 1);
    pct75 = prctile(tfp_mat, 75, 1);

    subplot(2, 3, v);
    plot(time, -pct25 * 100, '--', 'Color', dark_blue, 'LineWidth', 1); hold on;
    plot(time, -pct75 * 100, '--', 'Color', dark_blue, 'LineWidth', 1);
    plot(time, -mean_tfp * 100, '-',  'Color', dark_blue, 'LineWidth', 1.6);

    title([titles{v}, ' Response']);
    xlabel('Quarters'); ylabel('Percent');
    grid on; box off;

    y_max = 1.1 * max(abs([-pct25, -pct75] * 100), [], 'all');
    axis([0 T -y_max y_max]);

    if v == 1
        legend('25th pct.', '75th pct.', 'Mean', 'Location', 'best', 'FontSize', 9);
    end
end

exportgraphics(fig_tfp, 'figures/IRF_TFP_mean_percentiles.pdf', 'ContentType', 'vector');

% === Financial Shock: Mean + Percentile IRFs ===

fig_fin = figure('Name','Financial Shock: Mean + Percentiles', ...
                 'Color','w', 'Position', [100 100 1000 600]);

for v = 1:length(var_base)
    fin_mat = nan(N, T);
    varname_base = var_base{v};

    for i = 1:N
        fin_irf = sprintf('%s%d_eps_xi%d', varname_base, i, i);
        if isfield(irfs_het_xi, fin_irf)
            fin_mat(i, :) = irfs_het_xi.(fin_irf)(1:T);
        end
    end

    mean_fin = mean(fin_mat, 1, 'omitnan');
    pct25 = prctile(fin_mat, 25, 1);
    pct75 = prctile(fin_mat, 75, 1);

    subplot(2, 3, v);
    plot(time, -pct25 * 100, '--', 'Color', dark_green, 'LineWidth', 1); hold on;
    plot(time, -pct75 * 100, '--', 'Color', dark_green, 'LineWidth', 1);
    plot(time, -mean_fin * 100, '-',  'Color', dark_green, 'LineWidth', 1.6);

    title([titles{v}, ' Response']);
    xlabel('Quarters'); ylabel('Percent');
    grid on; box off;

    y_max = 1.1 * max(abs([-pct25, -pct75] * 100), [], 'all');
    axis([0 T -y_max y_max]);

    if v == 1
        legend('25th pct.', '75th pct.', 'Mean', 'Location', 'best', 'FontSize', 9);
    end
end

exportgraphics(fig_fin, 'figures/IRF_fin_mean_percentiles.pdf', 'ContentType', 'vector');


% Distributional Aggregate Response Visualization

%{
% Extract firm-level IRFs from Dynare output
N=250;
irf_length = 30;
timeRange = 0:irf_length-1; 
firm_vars = {'k', 'n', 'b', 'd', 'y'}; % Firm-specific variables
shocks = {'eps_z', 'eps_xi'}; % Aggregate shocks

% Initialize storage for computed moments
irf_means = struct();
irf_stds = struct();
irf_quantiles = struct();

% Compute statistics for each variable
for v = 1:length(firm_vars)
    var_name = firm_vars{v};
    
    for s = 1:length(shocks)
        shock_name = shocks{s};
        
        % Collect firm responses over time
        firm_responses = zeros(N, irf_length);
        for i = 1:N
            firm_responses(i, :) = oo_.irfs.([var_name num2str(i) '_' shock_name]) * 100;
        end
        
        % Compute moments
        irf_means.(var_name).(shock_name) = mean(firm_responses, 1);
        irf_stds.(var_name).(shock_name) = std(firm_responses, 1);
        irf_quantiles.(var_name).(shock_name) = prctile(firm_responses, [10, 25, 50, 75, 90], 1);
    end
end

% === Plot Improved Aggregate and Distributional IRFs ===
figure('Name', 'Aggregate and Distributional IRFs', ...
       'NumberTitle', 'off', ...
       'Position', [100, 100, 1400, 900], ...
       'Color', 'w');


firm_vars = {'k', 'n', 'b', 'd', 'y'};
shocks = {'eps_z', 'eps_xi'};

tiledlayout(length(firm_vars), length(shocks), ...
            'Padding', 'compact', 'TileSpacing', 'compact');

for v = 1:length(firm_vars)
    var_name = firm_vars{v};

    for s = 1:length(shocks)
        shock_name = shocks{s};
        nexttile;

        % Shaded 10–90 percentile
        fill([timeRange, fliplr(timeRange)], ...
             [irf_quantiles.(var_name).(shock_name)(1,:), ...
              fliplr(irf_quantiles.(var_name).(shock_name)(5,:))], ...
             [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
        hold on;

        % Shaded 25–75 percentile
        fill([timeRange, fliplr(timeRange)], ...
             [irf_quantiles.(var_name).(shock_name)(2,:), ...
              fliplr(irf_quantiles.(var_name).(shock_name)(4,:))], ...
             [0.2 0.5 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4); 

        % Median (dashed black)
        plot(timeRange, irf_quantiles.(var_name).(shock_name)(3,:), ...
             'k--', 'LineWidth', 1.25);

        % Mean (solid dark blue)
        plot(timeRange, irf_means.(var_name).(shock_name), ...
             'Color', [0 0.2 0.5], 'LineWidth', 1.75);

        % Aesthetics
        title([upper(var_name), ' response to ', shock_name], 'FontSize', 11);
        xlabel('Periods', 'FontSize', 9);
        ylabel('% Deviation', 'FontSize', 9);
        grid on;
        xlim([timeRange(1), timeRange(end)]);

        if v == 1 && s == 1
            legend({'10–90th', '25–75th', 'Median', 'Mean'}, ...
                   'FontSize', 8, 'Location', 'southoutside', ...
                   'Orientation', 'horizontal', 'Box', 'off');
        end
    end
end

%}