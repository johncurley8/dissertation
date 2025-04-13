function [ys,params,check] = het_xi_steadystate(ys,exo,M_,options_)

NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

check = 0;

N = 250; % Number of firms

rng(123);
sigma_xi = 0.2;
mu_xi = -0.5 * sigma_xi^2;
xi_vec = lognrnd(mu_xi, sigma_xi, [N,1]);
%xi_vec = ones(N, 1);
z = 1;

options = optimset('display','off');
xi_0 = 0.2338;

[xi, fval, exitflag] = fsolve(@(x)get_xi(N, x, xi_vec, z, tau, betta, delta, theta, nu, BY_ratio, ppsi), xi_0, options);
if exitflag < 1
    check = 1;
end

R = (1 - tau) / betta + tau;

mu_vec = (1 - R * betta) ./ (xi * xi_vec .* (R * (1 - tau) / (R - tau)));
%mu_vec = (1 - R * betta) * ones(N, 1) / (R * (1 - tau) / (R - tau));

% Define total labor
n_guess = 0.3;
n = n_guess;
% Initial guess for labor allocation across `N-1` firms (evenly distributed)
n_guess_vec = repmat(n / N, N-1, 1);

[n_solved, fval, exitflag] = fsolve(@(n_vec) compute_labor_residual(N, n_vec, z, xi, xi_vec, betta, delta, theta, nu, mu_vec, n), n_guess_vec, options);
if exitflag <= 0
    error('fsolve failed for labor allocation.');
end

n_firms = [n_solved; n - sum(n_solved)];
k_firms = zeros(N, 1);

for i = 1:N
    k_firms(i) = ((((1 - xi * xi_vec(i) * mu_vec(i)) / betta) - (1 - delta)) / ((1 - mu_vec(i)) * theta * z * n_firms(i)^nu))^(1 / (theta - 1));
end
k = sum(k_firms);

w = nu * sum(z .* k_firms.^theta .* n_firms.^(nu-1) .* (1 - mu_vec)) / N;

%b_firms = (ppsi .* z .* k_firms.^theta .* n_firms.^nu - (xi .* xi_vec) .* k_firms) * (tau - R) / (1 - tau);
b_firms = (ppsi .* z .* k_firms.^theta .* n_firms.^nu ./ (xi * xi_vec) - k_firms) * (tau - R) / (1 - tau);
b = sum(b_firms);

d_firms = (1 - delta) * k_firms + z .* k_firms.^theta .* n_firms.^nu - w .* n_firms - b_firms + b_firms / R - k_firms;
d = sum(d_firms);

c = w * n + b - b / R + d;
alppha = w / c^siggma * (1 - n);

y_firms = z .* k_firms.^theta .* n_firms.^nu;
y = sum(y_firms);

%% Additional variables
yhat = 0;
yhat1 = 0;
chat = 0;
v = d / (1 - betta);
byhat = 0;
dyhat = d / y;
vyhat = 0;
invest = delta * k;
r = (R - tau) / (1 - tau) - 1;
ihat = 0;
nhat = 0;

params = NaN(NumberOfParameters,1);
for iter = 1:length(M_.params)
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ]);
end

for i = 1:N
    eval(['xi', num2str(i), ' = xi_vec(i);']);
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
save('steady_state_firmdata_xi.mat', 'xi_vec', 'k_firms', 'n_firms', 'y_firms', ...
     'b_firms', 'd_firms', 'mu_vec');

end

function resid = get_xi(N, xi, xi_vec, z, tau, betta, delta, theta, nu, BY_ratio, ppsi)

    % Compute the steady-state return on capital
    R = (1 - tau) / betta + tau;
    r = (R - tau) / (1 - tau) - 1;
    
mu_vec = (1 - R * betta) ./ (xi * xi_vec .* (R * (1 - tau) / (R - tau)));
%mu_vec = (1 - R * betta) * ones(N, 1) / (R * (1 - tau) / (R - tau));

n_guess = 0.3;
n = n_guess;
n_guess_vec = repmat(n / length(xi_vec), length(xi_vec)-1, 1);

[n_solved, ~, exitflag] = fsolve(@(n_vec) compute_labor_residual(N, n_vec, z, xi, xi_vec, betta, delta, theta, nu, mu_vec, n), n_guess_vec, optimset('display','off'));
if exitflag <= 0
    error('fsolve failed for labor allocation.');
end

n_firms = [n_solved; n - sum(n_solved)];
k_firms = ((((1 - xi .* xi_vec .* mu_vec) / betta) - (1 - delta)) ./ ((1 - mu_vec) .* theta .* z .* n_firms.^nu)).^(1 / (theta - 1));

y_firms = z .* k_firms.^theta .* n_firms.^nu;
y = sum(y_firms);

b_firms = (ppsi .* y_firms ./ (xi .* xi_vec) - k_firms) * (tau - R) / (1 - tau);
%b_firms = (ppsi .* y_firms - (xi .* xi_vec) .* k_firms) * (tau - R) / (1 - tau);
b = sum(b_firms);

resid = BY_ratio - (b / (1 + r) / y);
end

function res = compute_labor_residual(N, n_vec, z, xi, xi_vec, betta, delta, theta, nu, mu_vec, n)
n_vec = max(real(n_vec), 1e-6);
n_residual = max(real(n - sum(n_vec)), 1e-6);
n_firms = [n_vec(:); n_residual];

k_firms = ((((1 - xi .* xi_vec .* mu_vec) / betta) - (1 - delta)) ./ ((1 - mu_vec) .* theta .* z .* n_firms.^nu)).^(1 / (theta - 1));

res = zeros(N-1, 1);
for i = 1:N-1
    res(i) = nu * z * k_firms(i)^theta * n_firms(i)^(nu - 1) * (1 - mu_vec(i)) - ...
             nu * z * k_firms(N)^theta * n_firms(N)^(nu - 1) * (1 - mu_vec(N));
end
end