/*
 * John Curley
 */
        
@#define N = 2

var 
    % Aggregate Variables
    c ${c}$ (long_name='Consumption')
    n ${n}$ (long_name='Total Labor')
    w ${w}$ (long_name='Wage')
    k ${k}$ (long_name='Total Capital')
    R ${R}$ (long_name='Effective Gross Interest Rate')
    d ${d}$ (long_name='Total Dividend')
    b ${b}$ (long_name='Total Bonds')

    % Firm-Specific Variables
@# for i in 1:N
    n@{i} ${n_@{i}}$ (long_name='Labor Type @{i}')
    k@{i} ${k_@{i}}$ (long_name='Capital Type @{i}')
    b@{i} ${b_@{i}}$ (long_name='Debt Type @{i}')
    d@{i} ${d_@{i}}$ (long_name='Dividend Type @{i}')
    y@{i} ${y_@{i}}$ (long_name='Output Type @{i}')
    mu@{i} ${\mu_@{i}}$ (long_name='Lagrange Multiplier Type @{i}')
    z@{i} ${z_@{i}}$ (long_name='Technology Type @{i}')
@# endfor

    % Quasi-exogenous
     z ${z}$ (long_name='Aggregate Technology ')
     xi ${\xi}$ (long_name='Financial Conditions')
    
    y ${y}$ (long_name='Output')
    yhat ${\hat y}$ (long_name='output deviation from trend')
    chat ${\hat c}$ (long_name='consumption deviation from trend') 
    byhat ${\hat y}$ (long_name='debt-repurchase to GDP ratio deviation from trend') 
    dyhat ${\hat y}$ (long_name='dividend to GDP ratio deviation from trend')
    v ${v}$ (long_name='Value of the Firm')
    vyhat ${\hat y}$ (long_name='Equity value deviation from trend')
    r ${r}$ (long_name='Gross Interest Rate')
    invest ${i}$ (long_name='Investment')
    ihat ${\hat i}$ (long_name='investment deviation from trend')
    nhat ${\hat n}$ (long_name='hours deviation from trend')
    % muhat ${\hat \mu}$ (long_name='Lagrange multiplier from trend')
         ;

varexo eps_z ${\varepsilon_z}$ (long_name='Aggregate technology shock')
       eps_xi ${\varepsilon_{\xi}}$ (long_name='Aggregate financial Shock')
@# for i in 1:N
    eps_z@{i} ${\varepsilon_{z@{i}}}$ (long_name='Technology Type @{i} shock')
@# endfor
         ;

parameters theta ${\theta}$ (long_name='capital share')
        nu ${\nu}$ (long_name='labor share')
        betta ${\beta}$ (long_name='discount factor')
        alppha ${\alpha}$ (long_name='disutility from work')
        delta ${\delta}$ (long_name='depreciation')
        tau ${\tau}$ (long_name='tax wedge')
        kappa ${\kappa}$ (long_name='equity cost')
        siggma ${\sigma}$ (long_name='risk aversion')
        sigma_z ${\sigma_z}$ (long_name='std_z')
        @# for i in 1:N
            sigma_z@{i}  % Firm-specific technology shock volatility
        @# endfor
        sigma_xi ${\sigma_xi}$ (long_name='std_xi')
        BY_ratio ${(\bar b/(1+\bar r)/\bar Y}$ (long_name='Debt output ratio')
        ppsi ${\psi}$ (long_name='psi')
        ppi ${\pi}$ (long_name='debt/dividends ratio')

        A11 ${A_{11}}$ (long_name='A_11')
        A12 ${A_{12}}$ (long_name='A_12')
        A21 ${A_{21}}$ (long_name='A_21')
        A22 ${A_{22}}$ (long_name='A_22')
        ;
     
    siggma=1;
    theta = 0.25;
    nu = 0.6;
    betta = 0.9825;
    delta = 0.025;
    tau = 0.35;        
    BY_ratio=3.36;
    ppsi = 1.0;

    kappa = 0.146;       
    kappa_store=kappa;
    sigma_xi = 0.0098;  
    sigma_z = 0.0055; 
    @# for i in 1:N
        sigma_z@{i} = 0.055;  % Firm-specific TFP shocks
    @# endfor
    A11 = 0.9457;
    A12 = 0;
    A21 = 0;
    A22 = 0.9703;
    options_.TeX=1;

model;
[name='FOC labor, equation 1 on p. 4 Appendix']
w/c^siggma-alppha/(1-n)=0;

[name='Euler equation, equation 2 on p. 4 Appendix']
c^(-siggma)=betta*((R-tau)/(1-tau))*c(+1)^(-siggma);

[name='Budget constraint household, equation 3 on p. 4 Appendix']
w*n+b(-1)-b/R+d-c=0;

@# for i in 1:N
    [name='FOC labor input firm @{i}']
    (nu) * z * z@{i} * k@{i}(-1)^theta * n@{i}^(nu-1) = w * (1 / (1 - mu@{i} * ppsi * (1 + 2 * kappa * (d@{i} - steady_state(d@{i})))));
@# endfor

@# for i in 1:N
    [name='FOC capital firm @{i}']
    betta * (c / c(+1))^siggma * ((1 + 2*kappa*(d@{i} - steady_state(d@{i}))) / (1 + 2*kappa*(d@{i}(+1) - steady_state(d@{i})))) *
    (1 - delta + (1 - mu@{i}(+1) * ppsi * (1 + 2*kappa*(d@{i}(+1) - steady_state(d@{i})))) * theta * z(+1) * z@{i}(+1) * k@{i}^(theta-1) * n@{i}(+1)^(nu))
    + xi * mu@{i} * (1 + 2*kappa*(d@{i} - steady_state(d@{i}))) = 1;
@# endfor

@# for i in 1:N
    [name='FOC bonds firm @{i}']
    R * betta * (c / c(+1))^siggma * ((1 + 2*kappa*(d@{i} - steady_state(d@{i}))) / (1 + 2*kappa*(d@{i}(+1) - steady_state(d@{i})))) 
    + xi * mu@{i} * (1 + 2*kappa*(d@{i} - steady_state(d@{i}))) * (R * (1 - tau) / (R - tau)) = 1;
@# endfor

@# for i in 1:N
    [name='Budget constraint firm @{i}']
    (1 - delta) * k@{i}(-1) + z * z@{i} * k@{i}(-1)^theta * n@{i}^nu - w * n@{i} - b@{i}(-1) + b@{i} / R - k@{i} - (d@{i} + kappa * (d@{i} - steady_state(d@{i}))^2) = 0;
@# endfor

@# for i in 1:N
    [name='Idiosyncratic TFP firm @{i}']
    log(z@{i} / steady_state(z@{i})) = 0.75 * log(z@{i}(-1) / steady_state(z@{i})) + A12 * log(xi(-1) / steady_state(xi)) + eps_z@{i};
@# endfor

@# for i in 1:N
    [name='Production function firm @{i}']
    y@{i} = z * z@{i} * k@{i}(-1)^theta * n@{i}^nu;
@# endfor

@# for i in 1:N
    [name='Enforcement constraint firm @{i}', bind = 'MU@{i}']
    xi * (k@{i} - b@{i} * ((1 - tau) / (R - tau))) = ppsi * z * z@{i} * k@{i}(-1)^theta * n@{i}^nu;
@# endfor

@# for i in 1:N
    [name='Enforcement constraint firm @{i}', relax = 'MU@{i}']
    mu@{i} = 0;
@# endfor

[name='Aggregate TFP process']
log(z/steady_state(z))=A11*log(z(-1)/steady_state(z))+eps_z; %+A12*log(xi(-1)/steady_state(xi))

[name='Aggregate financial conditions process']
log(xi/steady_state(xi))=A22*log(xi(-1)/steady_state(xi))+eps_xi;


%%%% MARKET CLEARING CONDITIONS (APPENDED BY mkt_clearing.m) %%%%

[name='Capital market clearing']
k = k1 + k2;

[name='Labor market clearing']
n = n1 + n2;

[name='Bond market clearing']
b = b1 + b2;

[name='Dividend market clearing']
d = d1 + d2;

[name='Output market clearing']
y = y1 + y2;

/*
end market clearing
*/

[name='Law of Motion of Capital']
invest=k-(1-delta)*k(-1);

[name='Definition output deviations from trend']
yhat=log(y)-log(steady_state(y));

[name='Definition consumption deviation from trend']
chat=log(c)-log(steady_state(c));

[name='Definition debt repurchase share in GDP']
byhat=(b(-1)/(1+r(-1))-b/(1+r))/y;    

[name='Definition equity payout to GDP ratio']
dyhat=d/y;

[name='Definition firm value']
v=d+betta*c/c(+1)*v(+1);

[name='Definition equity share']
vyhat=log(v/(k(-1)-b(-1)))-log(steady_state(v/(k(-1)-b(-1))));

[name='Definition before tax interest rate']
r=(R-tau)/(1-tau)-1;

[name='Definition investment deviation from trend']
ihat=log(invest)-log(steady_state(invest));

[name='Definition hours deviation from trend']
nhat=log(n)-log(steady_state(n));

end;

occbin_constraints;
@# for i in 1:N
    name 'MU@{i}';
    bind xi * (k@{i} - b@{i} * ((1 - tau) / (R - tau))) >= 1e-5 + ppsi * y@{i};
    relax xi * (k@{i} - b@{i} * ((1 - tau) / (R - tau))) + 1e-5 < ppsi * y@{i} - 1e-5;
@# endfor
end;

shocks;
    var eps_z= sigma_z^2;
    var eps_xi= sigma_xi^2;
@# for i in 1:N
    var eps_z@{i} = sigma_z@{i}^2;
@# endfor 
end;

write_latex_original_model;
write_latex_dynamic_model;
write_latex_parameter_table;

%% make VAR diagonal for generation of IRFS
        
set_param_value('A12',0);
set_param_value('A21',0);

%steady;


shocks (surprise);
var eps_z;
periods 1;
values .055; 
end;

occbin_setup(simul_periods = 30, simul_maxit = 1500,
simul_check_ahead_periods = 0,
simul_debug);

occbin_solver; // Solve the model with regime switching
occbin_graph(noconstant);

occbin_irfs.piecewise = oo_.occbin.simul.piecewise;
occbin_irfs.linear = oo_.occbin.simul.linear;
occbin_irfs.variable_names = M_.endo_names;

% Save the IRFs to a .mat file
save('occbin_irfs.mat', 'occbin_irfs');

/*
%% Display some statistics
fprintf('(b/(1+r)/Y: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))/(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact'))))/oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))
fprintf('Debt-Capital Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))/(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact'))))/oo_.dr.ys(strmatch('k',M_.endo_names,'exact')))
fprintf('Total Debt-Output Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))+oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))/(oo_.dr.ys(strmatch('y',M_.endo_names,'exact'))*(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact')))))
fprintf('Total Debt-Capital Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))+oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))/(oo_.dr.ys(strmatch('k',M_.endo_names,'exact'))*(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact')))))
fprintf('Capital-Output Ratio: %4.3f\n',(oo_.dr.ys(strmatch('k',M_.endo_names,'exact'))/oo_.dr.ys(strmatch('y',M_.endo_names,'exact'))))
*/

%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timeRange = 0:30;
/*

figure(1);
set(gcf, 'Name', 'Figure 1: IRF Comparison (First vs. Second Order)', 'NumberTitle', 'off');

subplot(2,4,1)
plot(-irfs_order1.yhat_eps_z * 100, 'b', 'LineWidth', 1.5) % First-order
hold on
plot(-irfs_order2.yhat_eps_z * 100, 'r--', 'LineWidth', 1.5) % Second-order
title('Output')

subplot(2,4,2)
plot(-irfs_order1.chat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.chat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Consumption')

subplot(2,4,3)
plot(-irfs_order1.nhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.nhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Hours')

subplot(2,4,4)
plot(-irfs_order1.ihat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.ihat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Investment')

subplot(2,4,5)
plot(-irfs_order1.byhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.byhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Debt rep./Y')

subplot(2,4,6)
plot(-irfs_order1.dyhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.dyhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Equity Payout/Y')

subplot(2,4,7)
plot(-irfs_order1.vyhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.vyhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Equity Value')

legend('First Order', 'Second Order')
%hold off


*/
/*
figure(1);
set(gcf, 'Name', 'Figure 1: Aggregate Impulse Responses', 'NumberTitle', 'off');
        subplot(2,4,1)
        plot(-oo_.irfs.yhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.yhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -0.8 0]);
        ll = legend('TFP shock', 'Financial Shock', 'Location', 'southwest');
        title('Output')
        
        subplot(2,4,2)
        plot(-oo_.irfs.chat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.chat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -0.6 0.06]);
        title('Consumption')
        
        subplot(2,4,3)
        plot(-oo_.irfs.nhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.nhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -1.2 0.8]);
        title('Hours')
        
        subplot(2,4,4)
        plot(-oo_.irfs.ihat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.ihat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -3.5 0.5]);
        title('Investment')
  
        subplot(2,4,5)
        plot(-oo_.irfs.byhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.byhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -1.9 2.4]);
        title('Debt rep./Y')
        
        subplot(2,4,6)
        plot(-oo_.irfs.dyhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.dyhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -1.5 1.5]);
        title('Equity Payout/Y')
        
        subplot(2,4,7)
        plot(-oo_.irfs.vyhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.vyhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -0.3 0.4]);
        title('Equity Value') 
*/
/*        
        subplot(2,4,8)
        plot(-oo_.irfs.muhat_eps_z*100,'b');
        hold on;
        plot(-oo_.irfs.muhat_eps_xi*100,'r--');
        title('Enforcement Constraint');
        legend('mu1, TFP1', 'mu1, Financial', 'Location', 'Best');
        
        subplot(2,4,8)
        plot(-oo_.irfs.muhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.muhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -20 25]);
        title('Multiplier, \mu')
*/
%filename = './Jermann_Quadrini_2012_RBC/graphs/Figure6_IRFs_alt';

% Save the figure as EPS (Encapsulated PostScript)
%print('-depsc2', filename);

%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
/*
figure(2);
set(gcf, 'Name', 'Figure 2: Firm-Level Responses to Aggregate TFP Shock', 'NumberTitle', 'off');
%clf;

% Capital Responses
subplot(2,3,1);
plot(-oo_.irfs.k1_eps_z*100, 'b'); % Firm 1 capital response to its TFP shock
hold on;
plot(-oo_.irfs.k2_eps_z*100, 'b-.'); % Firm 2 capital response to Firm 1's TFP shock
title('Capital (TFP Shock)');
legend('k1, TFP1', 'k2, TFP1', 'Location', 'Best');

% Labor Responses
subplot(2,3,2);
plot(-oo_.irfs.n1_eps_z*100, 'b');
hold on;
plot(-oo_.irfs.n2_eps_z*100, 'b-.');
title('Labor (TFP Shock)');
legend('n1, TFP1', 'n2, TFP1', 'Location', 'Best');

% Debt Responses
subplot(2,3,3);
plot(-oo_.irfs.b1_eps_z*100, 'b');
hold on;
plot(-oo_.irfs.b2_eps_z*100, 'b-.');
title('Debt (TFP Shock)');
legend('b1, TFP1', 'b2, TFP1', 'Location', 'Best');

% Dividend Responses
subplot(2,3,4);
plot(-oo_.irfs.d1_eps_z*100, 'b');
hold on;
plot(-oo_.irfs.d2_eps_z*100, 'b-.');
title('Dividends (TFP Shock)');
legend('d1, TFP1', 'd2, TFP1', 'Location', 'Best');

% Wage Responses
subplot(2,3,5);
plot(-oo_.irfs.w_eps_z*100, 'b');
title('Wage (TFP Shock)');
legend('w, TFP1', 'Location', 'Best');

% Firm-Specific Output Responses
subplot(2,3,6);
plot(-oo_.irfs.y1_eps_z*100, 'b');
hold on;
plot(-oo_.irfs.y2_eps_z*100, 'b-.');
title('Output (TFP Shock)');
legend('y1, TFP1', 'y2, TFP1', 'Location', 'Best');

for i = 1:6
    subplot(2,3,i)
    axis([0, 30, -2, 2]) % Adjust as needed based on IRF magnitude
end


%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
set(gcf, 'Name', 'Figure 3: Firm-Level Responses to Idiosyncratic TFP Shock', 'NumberTitle', 'off');
%clf;

% Capital Responses
subplot(2,3,1);
plot(-oo_.irfs.k1_eps_z1*100, 'b'); % Firm 1 capital response to its TFP shock
hold on;
plot(-oo_.irfs.k2_eps_z1*100, 'b-.'); % Firm 2 capital response to Firm 1's TFP shock
title('Capital (TFP Shock)');
legend('k1, TFP1', 'k2, TFP1', 'Location', 'Best');

% Labor Responses
subplot(2,3,2);
plot(-oo_.irfs.n1_eps_z1*100, 'b');
hold on;
plot(-oo_.irfs.n2_eps_z1*100, 'b-.');
title('Labor (TFP Shock)');
legend('n1, TFP1', 'n2, TFP1', 'Location', 'Best');

% Debt Responses
subplot(2,3,3);
plot(-oo_.irfs.b1_eps_z1*100, 'b');
hold on;
plot(-oo_.irfs.b2_eps_z1*100, 'b-.');
title('Debt (TFP Shock)');
legend('b1, TFP1', 'b2, TFP1', 'Location', 'Best');

% Dividend Responses
subplot(2,3,4);
plot(-oo_.irfs.d1_eps_z1*100, 'b');
hold on;
plot(-oo_.irfs.d2_eps_z1*100, 'b-.');
title('Dividends (TFP Shock)');
legend('d1, TFP1', 'd2, TFP1', 'Location', 'Best');

% Wage Responses
subplot(2,3,5);
plot(-oo_.irfs.w_eps_z1*100, 'b');
title('Wage (TFP Shock)');
legend('w, TFP1', 'Location', 'Best');

% Firm-Specific Output Responses
subplot(2,3,6);
plot(-oo_.irfs.y1_eps_z1*100, 'b');
hold on;
plot(-oo_.irfs.y2_eps_z1*100, 'b-.');
title('Output (TFP Shock)');
legend('y1, TFP1', 'y2, TFP1', 'Location', 'Best');

for i = 1:6
    subplot(2,3,i)
    axis([0, 30, -2, 2]) % Adjust as needed based on IRF magnitude
end

%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
set(gcf, 'Name', 'Figure 4: Firm-Level Responses to Financial Shock', 'NumberTitle', 'off');

% Capital Responses
subplot(2,3,1);
plot(-oo_.irfs.k1_eps_xi*100, 'r--'); % Firm 1 capital response to financial shock
hold on;
plot(-oo_.irfs.k2_eps_xi*100, 'r:'); % Firm 2 capital response to financial shock
title('Capital (Financial Shock)');
legend('k1, Financial', 'k2, Financial', 'Location', 'Best');

% Labor Responses
subplot(2,3,2);
plot(-oo_.irfs.n1_eps_xi*100, 'r--');
hold on;
plot(-oo_.irfs.n2_eps_xi*100, 'r:');
title('Labor (Financial Shock)');
legend('n1, Financial', 'n2, Financial', 'Location', 'Best');

% Debt Responses
subplot(2,3,3);
plot(-oo_.irfs.b1_eps_xi*100, 'r--');
hold on;
plot(-oo_.irfs.b2_eps_xi*100, 'r:');
title('Debt (Financial Shock)');
legend('b1, Financial', 'b2, Financial', 'Location', 'Best');

% Dividend Responses
subplot(2,3,4);
plot(-oo_.irfs.d1_eps_xi*100, 'r--');
hold on;
plot(-oo_.irfs.d2_eps_xi*100, 'r:');
title('Dividends (Financial Shock)');
legend('d1, Financial', 'd2, Financial', 'Location', 'Best');

% Wage Responses
subplot(2,3,5);
plot(-oo_.irfs.w_eps_xi*100, 'r--');
title('Wage (Financial Shock)');
legend('w, Financial', 'Location', 'Best');

% Firm-Specific Output Responses
subplot(2,3,6);
plot(-oo_.irfs.y1_eps_xi*100, 'r--');
hold on;
plot(-oo_.irfs.y2_eps_xi*100, 'r:');
title('Output (Financial Shock)');
legend('y1, Financial', 'y2, Financial', 'Location', 'Best');
*/


%%%%%%%%%%%%%%%%%%%%%%%% Do simulations for counterfactuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*

innovations=load('innovations_replication.mat','resid');           
emp_data=load('emp_data.mat');                                 % storing empirical data 

%% empirical data 
emp_GDP=emp_data.log_Real_GDP_detrended;    
emp_value_added=emp_data.log_Real_Business_value_added_detrended(:,1);
emp_equity=emp_data.equity_payout_detrended;
emp_debt=emp_data.debt_repurchases_detrended;
emp_hours=emp_data.log_Total_Private_hours_detrended;
emp_consumption=emp_data.log_Real_Personal_Consumption;
emp_investment=emp_data.log_Real_Investment_detrended;
timeline=emp_data.data_timeline;  
n_points=length(emp_GDP);      

%% Do "no frictions"-case first
        
set_param_value('kappa',0);
set_param_value('tau',0.00001);        % trick to leave mu determined in steady state (if set to zero, log of steady state of muhat (log(0)) is NaN)

    set_param_value('A12',-0.0091);
    set_param_value('A21',0.0321);
                               
stoch_simul(order=1,irf=105,nomoments,nocorr,nofunctions);

%% initialize IRF generation
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag);
shock_matrix = zeros(n_points,M_.exo_nbr);      %create shock matrix with number of time periods in colums

%% set shocks for 'productivity shocks only' 
        
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = zeros(1,n_points); 
y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_prod_nofric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);     % deviation from steady state  
        
%% set shocks for 'financial shocks only' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = zeros(1,n_points);
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_fin_nofric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);         
        
%% set shocks for 'both shocks' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_both_nofric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);


%% Now do model with frictions (baseline model)

    set_param_value('kappa',kappa_store);
    set_param_value('tau',0.35);

stoch_simul(order=1,bandpass_filter=[6,32],irf=105,nograph,nomoments,nocorr,nofunctions);

%% initialize IRF generation
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag);
shock_matrix = zeros(n_points,M_.exo_nbr);                
        
%% set shocks for 'productivity shocks only' 

shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = zeros(1,n_points); 
y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_prod_fric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);               
        
%% set shocks for 'technology shocks only' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = zeros(1,n_points);
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_fin_fric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);                    
        
%% set shocks for 'both shocks' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_both_fric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);     
        
*/
%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
figure('Name','Figure 2: Multiplier')

plot(timeline,100*y_IRF_both_fric(strmatch('muhat',M_.endo_names,'exact'),:));
axis([-inf inf -90 120]);
title('Lagrange Multiplier');
*/
%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
ff=figure('Name','Figure 3: Response to productivity shock only');

subplot(2,2,1)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('yhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*y_IRF_prod_nofric(strmatch('yhat',M_.endo_names,'exact'),:),'r--');
hh3=plot(timeline,100*emp_GDP,'g');
hold off
axis([-inf inf -14 8]);
plot_NBER_recessions([hh1;hh2;hh3]);     %plot shaded recession dates behind plot with handle h
ll=legend([hh1 hh2 hh3],'Baseline','No fin. fric.','Data');
set(ll,'Location','NorthWest');
title('GDP');

subplot(2,2,2)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('nhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*y_IRF_prod_nofric(strmatch('nhat',M_.endo_names,'exact'),:),'r--');
hh3=plot(timeline,100*emp_hours,'g');
hold off
axis([-inf inf -14 8]);        
plot_NBER_recessions([hh1;hh2;hh3]);     %plot shaded recession dates behind plot with handle h
title('Hours');

subplot(2,2,3)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('byhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_debt,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Debt repurchase');

subplot(2,2,4)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('dyhat',M_.endo_names,'exact'),:),'b-.');   
hold on
hh2=plot(timeline,100*emp_equity,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Equity payout');

%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff=figure('Name','Figure 4: Response to financial shocks only');

subplot(2,2,1)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('yhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_GDP,'g');
hold off
axis([-inf inf -14 8]);
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
ll=legend([hh1 hh2],'Model','GDP');
set(ll,'Location','NorthWest');
title('GDP');

subplot(2,2,2)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('nhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_hours,'g');
hold off
axis([-inf inf -14 8]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Hours');

subplot(2,2,3)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('byhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_debt,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Debt repurchase');

subplot(2,2,4)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('dyhat',M_.endo_names,'exact'),:),'b-.'); 
hold on
hh2=plot(timeline,100*emp_equity,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Equity payout');
        
equity_payout_filtered_model=bpf(100*y_IRF_both_fric(strmatch('dyhat',M_.endo_names,'exact'),:)',6,32,12);
fprintf('Std(Equity payout Model): %4.3f \n',nanstd(equity_payout_filtered_model))
equity_payout_filtered_data=bpf(100*emp_equity,6,32,12);
fprintf('Std(Equity payout Data): %4.3f \n',nanstd(equity_payout_filtered_data))


%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff=figure('Name','Figure 5: Response to both shocks');

subplot(2,2,1)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('yhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_GDP,'g');
hold off
ll=legend('Model','GDP');
set(ll,'Location','NorthWest');
axis([-inf inf -14 8]);
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('GDP');

subplot(2,2,2)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('nhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_hours,'g');
hold off
axis([-inf inf -14 8]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Hours');

subplot(2,2,3)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('byhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_debt,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Debt repurchase');

subplot(2,2,4)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('dyhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_equity,'g');
hold off
axis([-inf inf -12 15]);        
title('Equity payout');
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h

fprintf('Std(Y): \t %3.2f \t %3.2f \t %3.2f \t %3.2f\n',100*std(emp_GDP), 100*std(y_IRF_prod_fric(strmatch('yhat',M_.endo_names,'exact'),:)),100*std(y_IRF_fin_fric(strmatch('yhat',M_.endo_names,'exact'),:)),100*std(y_IRF_both_fric(strmatch('yhat',M_.endo_names,'exact'),:)))
fprintf('Std(N): \t %3.2f \t %3.2f \t %3.2f \t %3.2f\n',100*std(emp_hours), 100*std(y_IRF_prod_fric(strmatch('nhat',M_.endo_names,'exact'),:)),100*std(y_IRF_fin_fric(strmatch('nhat',M_.endo_names,'exact'),:)),100*std(y_IRF_both_fric(strmatch('nhat',M_.endo_names,'exact'),:)))
*/
