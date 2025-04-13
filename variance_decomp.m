% === Load saved data
load('var_decomp_xi.mat', 'var_decomp', 'shock_names', 'var_names');
load('var_decomp_tfp.mat', 'var_decomp', 'shock_names', 'var_names');

% === Check and fix if labels are flipped
[n_shocks, n_vars, ~] = size(var_decomp);

if length(shock_names) ~= n_shocks || length(var_names) ~= n_vars
    if length(shock_names) == n_vars && length(var_names) == n_shocks
        % Swap the labels
        temp = shock_names;
        shock_names = var_names;
        var_names = temp;
    else
        error('Mismatch between matrix dimensions and label sizes. Check var_decomp.mat contents.');
    end
end

% === Select horizon to display
h = 1;
slice = var_decomp(:, :, h);

% === Create labeled table
safe_var_names = matlab.lang.makeValidName(var_names);  % Sanitize for table
T = array2table(slice, ...
    'RowNames', shock_names, ...
    'VariableNames', safe_var_names);

% === Display
fprintf('Conditional variance decomposition at horizon %d:\n', h);
disp(T);
