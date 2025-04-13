% === Multi-Panel IRF Comparison: yhat, nhat, chat ===
load('occbin_irfs.mat');  % Load saved IRFs from OccBin

vars_to_plot = {'yhat', 'nhat', 'chat'};
titles = {'Output', 'Hours', 'Consumption'};
line_styles = {'b-', 'r--'};  % Linear (blue solid), Piecewise (red dashed)

figure('Name', 'Linear vs. OccBin IRFs', 'Color', 'w', 'Position', [100 100 900 300]);

for i = 1:length(vars_to_plot)
    varname = vars_to_plot{i};
    idx = strcmp(occbin_irfs.variable_names, varname);

    linear_response = occbin_irfs.linear(:, idx);
    piecewise_response = occbin_irfs.piecewise(:, idx);
    time = 1:length(linear_response);

    subplot(1,3,i);
    plot(time, -linear_response, line_styles{1}, 'LineWidth', 1.5); hold on;
    plot(time, -piecewise_response, line_styles{2}, 'LineWidth', 1.5);
    title(titles{i});
    xlabel('Quarters');
    ylabel('Percent');

    grid on; box off;
    axis([0 max(time) -0.01 0.01]);

    % Force scientific notation with exponent at top of axis
    ax = gca;
    ax.YRuler.Exponent = -2;
    ax.YAxis.ExponentMode = 'manual';

    if i == 1
        legend('Linear', 'Piecewise Linear', 'Location', 'northwest', 'FontSize', 9);
    end
end

exportgraphics(gcf, '../figures/IRF_occbin_compare.pdf', 'ContentType', 'vector');
