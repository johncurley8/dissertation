function [ys,params,check] = het_TFP_steadystate(ys,exo,M_,options_)

NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

check = 0;

N = 250; % Number of firms

% Set steadystate TFP
z = 1;

rng(123);
%rng(999982);
mu_log_z = -0.5 * sigma_log_z^2;
z_vec = lognrnd(mu_log_z, sigma_log_z, [N, 1]);
%z_vec = ones(N, 1);
save('steady_state_firmdata.mat', 'z_vec');

options = optimset('Display','off','TolFun',1e-6,'TolX',1e-8,'MaxIter',1e4);

% Solve for xi that satisfies calibration targets
xi_0 = 0.40338;
[xi, fval, exitflag] = fsolve(@(x)get_xi(x, N, z, z_vec, tau, betta, delta, theta, nu, BY_ratio, ppsi), xi_0, options);
if exitflag < 1
    check = 1;
end

R = (1 - tau) / betta + tau;

mu_vec = (1 - R * betta) ./ (xi * (R * (1 - tau) / (R - tau)) * ones(N, 1));
%mu_vec = (1 - R * betta) ./ ((R * (1 - tau) / (R - tau)) * ones(N, 1));

% Define total labor
n_guess = 0.3;
n = n_guess;

% Initial guess for labor allocation
n_guess_vec = z_vec(1:N-1) / sum(z_vec) * n;  % proportional to TFP

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
%b_firms = (ppsi .* z .* z_vec .* k_firms.^theta .* n_firms.^nu - xi * k_firms) * (tau - R) / (1 - tau);
b = sum(b_firms);


% Compute firm-specific dividends
d_firms = (1 - delta) * k_firms + z .* z_vec .* k_firms.^theta .* n_firms.^nu ...
          - w .* n_firms - b_firms + b_firms / R - k_firms;
d = sum(d_firms);

% Consumption
c = w * n + b - b / R + d;
alppha = w / c^siggma * (1 - n);


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

%% End model equations


params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params)
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
    eval(sprintf('byhat%d = 0;', i));
    eval(sprintf('dyhat%d = d_firms(%d) / y_firms(%d);', i, i, i));
    eval(sprintf('byhat%d = 0;', i));
    eval(sprintf('invest%d = delta * k_firms(%d);', i, i));
    eval(sprintf('ihat%d = 0;', i));
    eval(sprintf('nhat%d = 0;', i));
    eval(sprintf('yhat%d = 0;', i));
end

NumberOfEndogenousVariables = M_.orig_endo_nbr;
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end

% Save firm-level steady state data for external analysis/plots
save('steady_state_firmdata_TFP.mat', 'z_vec', 'k_firms', 'n_firms', 'y_firms', ...
     'b_firms', 'd_firms', 'mu_vec');
end
function resid = get_xi(xi, N, z, z_vec, tau, betta, delta, theta, nu, BY_ratio, ppsi)
    
    % Compute the steady-state return on capital
    R = (1 - tau) / betta + tau;
    r = (R - tau) / (1 - tau) - 1;
    
    % Compute firm-specific financial friction multipliers
    mu_vec = (1 - R * betta) ./ (xi * (R * (1 - tau) / (R - tau)) * ones(N, 1));
    %mu_vec = (1 - R * betta) ./ ((R * (1 - tau) / (R - tau)) * ones(N, 1));
    n_guess = 0.3; 
    n = n_guess; 

    % Solve for firm-specific labor allocations
    options = optimset('Display','off','TolFun',1e-6,'TolX',1e-8,'MaxIter',1e4);

    n_guess_vec = z_vec(1:N-1) / sum(z_vec) * n;  % proportional to TFP

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
    %b_firms = (ppsi .* y_firms - xi .* k_firms) * (tau - R) / (1 - tau);
    b = sum(b_firms);

    % Compute residual (deviation of debt-to-GDP from target)
    resid = (BY_ratio - (b / (1 + r) / y));

    % Check borrowing constraints
    constraint_firms = xi * (k_firms - b_firms * (1 - tau) / (R - tau));

    for i = 1:N
        disp(['Enforcement Constraint Firm ', num2str(i), ': ', num2str(constraint_firms(i)), ' =? ', num2str(y_firms(i))]);
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

    MPn = @(z_i, k_i, n_i, mu_i) ...
    log(nu) + log(z) + log(z_i) + theta * log(k_i) + (nu - 1) * log(n_i) + log(1 - mu_i);

    % Compute residual equations for labor market (ensuring proper indexing)
    res = zeros(N-1, 1);
    for i = 1:N-1
       res(i) = MPn(z_vec(i), k_firms(i), n_firms(i), mu_vec(i)) - ...
         MPn(z_vec(N), k_firms(N), n_firms(N), mu_vec(N));
    end

end

