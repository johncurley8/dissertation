function [ys,params,check] = het_TFP_Occ_steadystate(ys,exo,M_,options_)

NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

check = 0;

N = 2; % Number of firms

% Set steadystate TFP
z = 1;

sigma_log_z = 0.3;
mu_log_z = -0.5 * sigma_log_z^2;
z_vec = lognrnd(mu_log_z, sigma_log_z, [N, 1]);
%z_vec = ones(N, 1);

binding_regime = 0; % 0 means nonbinding

if isfield(M_, 'occbin') && isfield(M_.occbin, 'binding_regime') && M_.occbin.binding_regime == 0
    binding_regime = 1; 
end
disp(['Binding regime detected: ', num2str(binding_regime)]);

if binding_regime  %Solve using get_xi 

options = optimset('display', 'off');

% Solve for xi that satisfies calibration targets
xi_0 = 0.40338;
[xi, fval, exitflag] = fsolve(@(x)get_xi(x, N, z, z_vec, tau, betta, delta, theta, nu, BY_ratio, ppsi), xi_0, options);
if exitflag < 1
    check = 1;
end

R = (1 - tau) / betta + tau;
%R = 1/betta;

mu_vec = (1 - R * betta) ./ (xi * (R * (1 - tau) / (R - tau)) * ones(N, 1));

% Define total labor
n_guess = 0.3;
n = n_guess;

% Initial guess for labor allocation across `N-1` firms (evenly distributed)
n_guess_vec = repmat(n / N, N-1, 1);  

% Solve for labor allocations numerically using vectorized function
[n_solved, fval, exitflag] = fsolve(@(n_vec) compute_labor_residual(N, n_vec, z, z_vec, xi, betta, delta, theta, nu, mu_vec, n, ppsi), ...
                                     n_guess_vec, options);

% Check if `fsolve` converged
if exitflag <= 0
    error('fsolve failed to find a solution for firm-specific labor allocation.');
end

% Store solutions: Last firm's labor inferred from market clearing
n_firms = [n_solved; n - sum(n_solved)];

for i = 1:N
    eval(['n', num2str(i), ' = n_firms(i);']);
end


k_firms = zeros(N, 1);
for i = 1:N
    k_firms(i) = ((((1 - xi * mu_vec(i)) / betta) - (1 - delta)) / ...
                 ((1 - mu_vec(i)) * theta * z * z_vec(i) * n_firms(i)^nu))^(1 / (theta - 1));
end
k = sum(k_firms);

% Compute wage
w = nu * sum(z * z_vec .* k_firms.^theta .* n_firms.^(nu-1) .* (1 - mu_vec)) / N;

% Compute firm-specific bonds
b_firms = (ppsi .* z .* z_vec .* k_firms.^theta .* n_firms.^nu ./ xi - k_firms) * (tau - R) / (1 - tau);
b = sum(b_firms);


% Compute firm-specific dividends
d_firms = (1 - delta) * k_firms + z .* z_vec .* k_firms.^theta .* n_firms.^nu ...
          - w .* n_firms - b_firms + b_firms / R - k_firms;
d = sum(d_firms);

% Consumption
c = w * n + b - b / R + d;
alppha = w / c^siggma * (1 - n);

% Output
y_firms = zeros(N, 1); 
for i = 1:N
    y_firms(i) = z * z_vec(i) * k_firms(i)^theta * n_firms(i)^nu;
end
y = sum(y_firms);

% Additional variables
yhat = 0;
chat = 0;
v = d / (1 - betta);
byhat = 0;
dyhat = d / y;
vyhat = 0;
invest = delta * k;
r = (R - tau) / (1 - tau) - 1;
ihat = 0;
nhat = 0;
%muhat = 0;

%% End model equations

    params=NaN(NumberOfParameters,1);
    for iter = 1:length(M_.params) %update parameters set in the file
      eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
    end

    % Translate variables from vector form
    for i = 1:N
        eval(['z', num2str(i), ' = z_vec(i);']); 
        eval(['mu', num2str(i), ' = mu_vec(i);']);
        eval(['k', num2str(i), ' = k_firms(i);']);
        eval(['b', num2str(i), ' = b_firms(i);']);
        eval(['d', num2str(i), ' = d_firms(i);']);
        eval(['n', num2str(i), ' = n_firms(i);']);
        eval(['y', num2str(i), ' = y_firms(i);']);
    end

    NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
    for ii = 1:NumberOfEndogenousVariables
      varname = M_.endo_names{ii};
      eval(['ys(' int2str(ii) ') = ' varname ';']);
    end
disp('-------------------------');
disp(['R: ', num2str(R)]);
disp(['mu1: ', num2str(mu1)]);
disp(['mu2: ', num2str(mu2)]);
disp(['k1: ', num2str(k1)]);
disp(['k2: ', num2str(k2)]);
disp(['b1: ', num2str(b1)]);
disp(['b2: ', num2str(b2)]);
disp('-------------------------');

else 
    disp('Skipping get_xi for NBR steady state...');
    N=2;
    mu_vec = zeros(N,1);
    tau = 0;
    xi = 0.27;
    % Compute the steady-state return on capital
    R = (1 - tau) / betta + tau;
    r = (R - tau) / (1 - tau) - 1;
    %R = 1/betta;

    % Define total labor
    n_guess = 0.3;
    n = n_guess;
    
    % Initial guess for labor allocation across `N-1` firms (evenly distributed)
    n_guess_vec = repmat(n / N, N-1, 1);  

    options = optimset('display', 'off');
    % Solve for labor allocations numerically using vectorized function
    [n_solved, fval, exitflag] = fsolve(@(n_vec) compute_labor_residual(N, n_vec, z, z_vec, xi, betta, delta, theta, nu, mu_vec, n, ppsi), ...
                                         n_guess_vec, options);
    
    % Check if `fsolve` converged
    if exitflag <= 0
        error('fsolve failed to find a solution for firm-specific labor allocation.');
    end
    
    % Store solutions: Last firm's labor inferred from market clearing
    n_firms = [n_solved; n - sum(n_solved)];
    
    for i = 1:N
        eval(['n', num2str(i), ' = n_firms(i);']);
    end
    
    k_firms = zeros(N, 1);
    for i = 1:N
        k_firms(i) = ((((1 - xi * mu_vec(i)) / betta) - (1 - delta)) / ...
                     ((1 - mu_vec(i)) * theta * z * z_vec(i) * n_firms(i)^nu))^(1 / (theta - 1));
    end
    k = sum(k_firms);
    
    % Compute wage
    w = nu * sum(z * z_vec .* k_firms.^theta .* n_firms.^(nu-1) .* (1 - mu_vec)) / N;
    disp(['Computed wage (w): ', num2str(w)]);

    % Solve for firm-level bonds and dividends
    [b_firms, d_firms] = solve_nonbinding_debt(N, betta, BY_ratio, delta, theta, nu, z, z_vec, n_firms, k_firms, w, R, r, ppsi, xi, tau);
    
    % Aggregate
    b = sum(b_firms);
    d = sum(d_firms);

    % Output
    y_firms = zeros(N, 1); 
    for i = 1:N
        y_firms(i) = z * z_vec(i) * k_firms(i)^theta * n_firms(i)^nu;
    end
    y = sum(y_firms);

    % Print individual components
    % disp(['Computed Firm 1 Debt (b1): ', num2str(b1)]);
    % disp(['Computed Firm 2 Debt (b2): ', num2str(b2)]);
    % disp(['Computed Firm 1 Dividends (d1): ', num2str(d1)]);
    % disp(['Computed Firm 2 Dividends (d2): ', num2str(d2)]);
    % disp(['Total Dividends (d1 + d2): ', num2str(d)]);


    % Compute consumption
    c = w * n + b - b / R + d;
    alppha = w / c^siggma * (1 - n);

    % Additional variables
    yhat = 0;
    chat = 0;
    v = d / (1 - betta);
    byhat = 0;
    dyhat = d / y;
    vyhat = 0;
    invest = delta * k;
    r = (R - tau) / (1 - tau) - 1;
    ihat = 0;
    nhat = 0;

    % Compute total household inflows and outflows
    hh_income = w * n + b - b / R + d;
    hh_expenditures = c;
    hh_budget_residual = hh_income - hh_expenditures;

    % Check if household budget constraint holds
    if abs(hh_budget_residual) > 1e-6
       error('Nonbinding regime: Household budget constraint mismatch detected.');
    end

    
    params=NaN(NumberOfParameters,1);
    for iter = 1:length(M_.params) 
        eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ]);
    end
    
    % Translate variables from vector form
    for i = 1:N
        eval(['z', num2str(i), ' = z_vec(i);']); 
        eval(['mu', num2str(i), ' = mu_vec(i);']);
        eval(['k', num2str(i), ' = k_firms(i);']);
        eval(['b', num2str(i), ' = b_firms(i);']);
        eval(['d', num2str(i), ' = d_firms(i);']);
        eval(['n', num2str(i), ' = n_firms(i);']);
        eval(['y', num2str(i), ' = y_firms(i);']);
    end

    NumberOfEndogenousVariables = M_.orig_endo_nbr; 
    for ii = 1:NumberOfEndogenousVariables
        varname = M_.endo_names{ii};
        eval(['ys(' int2str(ii) ') = ' varname ';']);
    end

    % Assign `ys` to Dynare's steady state vector
    oo_.steady_state = ys;

end
end

function [b_firms, d_firms] = solve_nonbinding_debt(N, betta, BY_ratio, delta, theta, nu, z, z_vec, n_firms, k_firms, w, R, r, ppsi, xi, tau)
    % Output
    y_firms = zeros(N, 1); 
    for i = 1:N
        y_firms(i) = z * z_vec(i) * k_firms(i)^theta * n_firms(i)^nu;
    end
    y = sum(y_firms);
    
    b = BY_ratio * y * (1+r);

    for i = 1:N
        % Debt share by output
        b_firms(i) = (y_firms(i) / y) * b;
    
        % Dividends
        d_firms(i) = -delta * k_firms(i) + y_firms(i) - w * n_firms(i) - b_firms(i) + b_firms(i) / R;
    
        % Binding debt level
        b_firms_bind(i) = (ppsi * y_firms(i) / xi - k_firms(i)) * (tau - R) / (1 - tau);
    end
    
    % Check that nonbinding debt is lower than binding debt
    for i = 1:N
        if b_firms(i) > b_firms_bind(i)
            warning(['Nonbinding regime debt (b' num2str(i) ') is higher than binding regime: ' ...
                     num2str(b_firms(i)) ' > ' num2str(b_firms_bind(i))]);
        end
    end

    % Debugging prints
    % disp(['Computed Firm 1 Debt (b1): ', num2str(b1), ' | Binding Regime (b1_bind): ', num2str(b1_bind)]);
    % disp(['Computed Firm 2 Debt (b2): ', num2str(b2), ' | Binding Regime (b2_bind): ', num2str(b2_bind)]);
    % disp(['Computed Firm 1 Dividends (d1): ', num2str(d1)]);
    % disp(['Computed Firm 2 Dividends (d2): ', num2str(d2)]);

    disp(['BY_ratio: ', num2str(BY_ratio)]);
    disp(['Output: ', num2str(y)]);
    disp(['Debt/Output ratio: ', num2str((b/(1+r))/y)]);
end

function resid = get_xi(xi, N, z, z_vec, tau, betta, delta, theta, nu, BY_ratio, ppsi)
    
    % Compute the steady-state return on capital
    R = (1 - tau) / betta + tau;
    r = (R - tau) / (1 - tau) - 1;
    
    % Compute firm-specific financial friction multipliers
    mu_vec = (1 - R * betta) ./ (xi * (R * (1 - tau) / (R - tau)) * ones(N, 1));

    % Define total labor (initial guess)
    n_guess = 0.3; 
    n = n_guess; 

    % Solve for firm-specific labor allocations
    options = optimset('Display', 'off');
    n_guess_vec = repmat(n / N, N-1, 1);

    [n_solved, fval, exitflag] = fsolve(@(n_vec) compute_labor_residual(N, n_vec, z, z_vec, xi, betta, delta, theta, nu, mu_vec, n, ppsi), ...
                                         n_guess_vec, options);

    % Store solutions
    n_firms = [n_solved; n - sum(n_solved)]; % Enforce labor market clearing

    % Solve for firm-specific capital
    k_firms = ((((1 - xi * mu_vec) / betta) - (1 - delta)) ./ ...
               ((1 - mu_vec) .* theta .* z .* z_vec .* n_firms.^nu)).^(1 / (theta - 1));

    % Compute firm-specific output
    y_firms = z .* z_vec .* k_firms.^theta .* n_firms.^nu;
    y = sum(y_firms);

    % Compute firm-specific bonds
    b_firms = (ppsi .* y_firms / xi - k_firms) * (tau - R) / (1 - tau);
    b = sum(b_firms);

    % Compute residual (deviation of debt-to-GDP from target)
    resid = (BY_ratio - (b / (1 + r) / y));

    % Check borrowing constraints
    constraint_firms = xi * (k_firms - b_firms * (1 - tau) / (R - tau));

    for i = 1:N
        disp(['Enforcement Constraint Firm ', num2str(i), ': ', num2str(constraint_firms(i)), ' =? ', num2str(ppsi * y_firms(i))]);
    end
end

function res = compute_labor_residual(N, n_vec, z, z_vec, xi, betta, delta, theta, nu, mu_vec, n, ppsi)
    global iter_count;
    iter_count = iter_count + 1;

    n_vec = max(real(n_vec), 1e-6);
    n_residual = max(real(n - sum(n_vec)), 1e-6);

    n_firms = [n_vec(:); n_residual];

    k_firms = zeros(N, 1);
    for i = 1:N
        numerator = ((1 - xi * mu_vec(i)) / betta) - (1 - delta);
        denominator = (1 - mu_vec(i)) * theta * z * z_vec(i) * (n_firms(i).^nu);  
        numerator = max(real(numerator), 1e-6); 
        denominator = max(real(denominator), 1e-6); 

        k_firms(i) = (numerator / denominator)^(1 / (theta - 1));
    end

    % Compute residual equations for labor market (ensuring proper indexing)
    res = zeros(N-1, 1);
    for i = 1:N-1
        res(i) = nu * z * z_vec(i) * k_firms(i)^theta * n_firms(i)^(nu - 1) * (1 - mu_vec(i)) - ...
                 nu * z * z_vec(N) * k_firms(N)^theta * n_firms(N)^(nu - 1) * (1 - mu_vec(N));
    end

end
