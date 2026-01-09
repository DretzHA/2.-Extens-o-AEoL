clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES GERAIS
%  ========================================================================
conf.k = 100;          
conf.delta = 0.2;    
conf.velocities = 1*ones(10, 1); 

conf.Symbol_Duration = 1; 

% mu = service rate | rho = system load | lambda = arrival rate
% Eq. 4 - Timely Updates over an Erasure Channel
avg_service_time = conf.k / (1 - conf.delta);
avg_service_time_sec = avg_service_time * conf.Symbol_Duration;

mu_eff = 1 / avg_service_time_sec; 

conf.rho_total_vec = 0.1:1:100; 
conf.lambda_total_vec = conf.rho_total_vec * mu_eff; 

conf.N_updates = 200;   
conf.N_MC = 3000;

user_scenarios = [1]; 

%% ========================================================================
%  2. Simulation
%  ========================================================================

for i = 1:length(user_scenarios)
    n_u = user_scenarios(i);
    conf.num_users = n_u;

    if length(conf.velocities) < n_u
        conf.velocities = ones(n_u, 1); 
    end
    v = conf.velocities(1);
    
    fprintf('Processing N=%d users...\n', n_u);
   
    % Monte Carlo
    [AEoL_NP_Sim, AEoL_P_Sim, AEoL_Sys_NP, AEoL_Sys_P] = sim_mc_multi_erasure(conf);
    
    % ---------------------------------------------------------------------
    % Theory
    % ---------------------------------------------------------------------
    aeol_theory_np = zeros(size(conf.lambda_total_vec));
    aeol_theory_p  = zeros(size(conf.lambda_total_vec));
    
    for j = 1:length(conf.lambda_total_vec)
        lambda_val = conf.lambda_total_vec(j); 
        lambda_per_slot = lambda_val * conf.Symbol_Duration;
        
        % --- NP Theory (Slots) ---
        num = lambda_per_slot * conf.k * (conf.k + conf.delta);
        denom = 2 * (1 - conf.delta) * (lambda_per_slot * conf.k + 1 - conf.delta);
        val_np_slots = (1/lambda_per_slot) + (conf.k/(1-conf.delta)) + (num/denom);
        aeol_theory_np(j) = v * val_np_slots * conf.Symbol_Duration;
        
        % --- P Theory (Slots) ---
        term_base = (exp(lambda_per_slot) - conf.delta) / (1 - conf.delta);
        val_p_slots = (1/lambda_per_slot) * (term_base ^ conf.k);
        aeol_theory_p(j) = v * val_p_slots * conf.Symbol_Duration;
    end

    %% ====================================================================
    %  Figures
    %  ===================================================================
    fig_name = sprintf('Combined Results (NP vs P) - %d Users', n_u);
    figure('Name', fig_name, 'Color', 'w', 'Position', [100 100 1200 700]);
    
    hold on; 
    grid on; grid minor; box on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.1, 'TickLabelInterpreter', 'latex');
    
    color_sys_p  = [0 0.4470 0.7410];      
    color_sys_np = [0.8500 0.3250 0.0980]; 
    color_theo   = [0.15 0.15 0.15];     
    
    user_colors = lines(n_u); 
    
    % Theory NP
    plot(conf.lambda_total_vec, aeol_theory_np, ...
        'LineStyle', '-', 'Color', color_theo, 'LineWidth', 2.5, ...
        'DisplayName', 'Theory (NP)');
        
    % Theory P
    plot(conf.lambda_total_vec, aeol_theory_p, ...
        'LineStyle', '--', 'Color', color_theo, 'LineWidth', 2.5, ...
        'DisplayName', 'Theory (P)');

    for u = 1:n_u
        col_u = user_colors(u, :);
        % User NP 
        plot(conf.lambda_total_vec, AEoL_NP_Sim(u, :), ...
            'LineStyle', '-', 'Color', [col_u 0.3], 'LineWidth', 1.0, ... 
            'HandleVisibility', 'off'); 
        % User P 
        plot(conf.lambda_total_vec, AEoL_P_Sim(u, :), ...
            'LineStyle', ':', 'Color', [col_u 0.6], 'LineWidth', 1.0, ...
            'HandleVisibility', 'off'); 
    end
    
    % Legend "Dummy plots" 
    plot(NaN, NaN, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'DisplayName', 'Indiv. Users (NP)');
    plot(NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'DisplayName', 'Indiv. Users (P)');

    % System NP
    plot(conf.lambda_total_vec, AEoL_Sys_NP, ...
        'LineStyle', '-', 'Color', color_sys_np, 'LineWidth', 2.2, ...
        'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', color_sys_np, ...
        'DisplayName', 'SYSTEM (NP)');

    % System P
    plot(conf.lambda_total_vec, AEoL_Sys_P, ...
        'LineStyle', '--', 'Color', color_sys_p, 'LineWidth', 2.2, ...
        'Marker', 's', 'MarkerSize', 6, 'MarkerFaceColor', color_sys_p, ...
        'DisplayName', 'SYSTEM (P)');

    xlabel('Total Update Rate $\lambda_{tot}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Average AEoL $\bar{\Delta}$', 'Interpreter', 'latex', 'FontSize', 14);
    
    title_str = sprintf('Scenario: %d Users (Shared Object) | k=%d, $\\delta=%.1f$', ...
                  n_u, conf.k, conf.delta);
    title(title_str, 'Interpreter', 'latex', 'FontSize', 14);
    
    lgd = legend('Location', 'northeast', 'NumColumns', 2, 'FontSize', 10);
    lgd.ItemTokenSize = [20, 18]; 
    
    yl = ylim;
    ylim([0 yl(2)*1.05]); 
    xlim([min(conf.lambda_total_vec) max(conf.lambda_total_vec)]);
    
    drawnow;
end