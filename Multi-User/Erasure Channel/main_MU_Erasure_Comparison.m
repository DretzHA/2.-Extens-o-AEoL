clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES GERAIS
%  ========================================================================
% Channel Parameters
conf.k = 16;          % Symbols per update
conf.delta = 0.1;    % Erasure Probability
conf.velocities = ones(10, 1); % Users Velocity
conf.Symbol_Duration = 1; % Ex: 1 ms para transmitir 1 símbolo (ou pacote base)


% mu = service rate | rho = system load | lambda = arrival rate
% Eq. 4 - Timely Updates over an Erasure Channel
avg_service_time = conf.k / (1 - conf.delta);
avg_service_time_sec = avg_service_time * conf.Symbol_Duration;

mu_eff = 1 / avg_service_time; % service rate

conf.rho_total_vec = 0.1:0.2:5; % system load
conf.lambda_total_vec = conf.rho_total_vec * mu_eff; % arrival rate

conf.N_updates = 100;   % Number of updates to be sent
conf.N_MC = 1000;        % Monte Carlo loops

user_scenarios = [1 2 3]; % Number of users in each scenario
colors = lines(length(user_scenarios) + 1); % color for figures


%% ========================================================================
%  2. Simulation
%  ========================================================================

for i = 1:length(user_scenarios)
    n_u = user_scenarios(i);
    conf.num_users = n_u;

    % Garante que tenhamos velocidades suficientes
    if length(conf.velocities) < n_u
        conf.velocities = ones(n_u, 1); 
    end
    current_velocities = conf.velocities(1:n_u);
    
    fprintf('Processing N=%d users...\n', n_u);
   

    [AEoL_NP_Sim, AEoL_P_Sim] = sim_mc_multi_erasure(conf);

    figure('Name', sprintf('Cenário %d Usuários', n_u), 'Color', 'w', 'Position', [100 100 1000 600]);
    hold on; grid on; box on;
    set(gca, 'FontSize', 11);
    
    user_colors = lines(n_u);

    for u = 1:n_u
        v_u = current_velocities(u);
        color_u = user_colors(u, :);
        
        % Non-Preemptive Simulation
        plot(conf.lambda_total_vec, AEoL_NP_Sim(u, :), ...
            'LineStyle', '-', 'Marker', 'x', 'Color', color_u, 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Sim NP U%d', u));

        % Preemptive Simulation
        plot(conf.lambda_total_vec, AEoL_P_Sim(u, :), ...
            'LineStyle', '--', 'Marker', 'o', 'Color', color_u, 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'DisplayName', sprintf('Sim P U%d', u));
  
    end

    aeol_theory_np = zeros(size(conf.lambda_total_vec));
    aeol_theory_p  = zeros(size(conf.lambda_total_vec));
    
    for j = 1:length(conf.lambda_total_vec)
        lambda_tot = conf.lambda_total_vec(j);
        lambda_u   = lambda_tot / n_u; 
        
        % Teoria Non-Preemptive
        % Eq (5)
        num = lambda_u * conf.k * (conf.k + conf.delta);
        denom = 2 * (1 - conf.delta) * (lambda_u * conf.k + 1 - conf.delta);
        val_np = (1/lambda_u) + (conf.k/(1-conf.delta)) + (num/denom);
        aeol_theory_np(j) = val_np * v_u; 
        
        % Teoria Preemptive
        % Eq (17) 
        term_base = (exp(lambda_u) - conf.delta) / (1 - conf.delta);
        val_p = (1/lambda_u) * (term_base ^ conf.k);
        aeol_theory_p(j) = val_p * v_u;
    end
        
    plot(conf.lambda_total_vec, aeol_theory_np, ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 2, ...
        'DisplayName', sprintf('Theo NP U%d', u));

    % Opcional: Plotar teoria P também, ou apenas NP para limpar o gráfico
    plot(conf.lambda_total_vec, aeol_theory_p, ...
       'LineStyle', '--', 'Color', 'k', 'LineWidth', 1, ...
       'DisplayName', sprintf('Theo P U%d', u));
    
    
    %% ====================================================================
    %  FORMATAÇÃO FINAL DA FIGURA
    %  ====================================================================
    xlabel('Total Update Rate (\lambda_{tot})');
    ylabel('AEoL (Average Age \times v)');
    title(sprintf('Cenário com %d Usuário(s) | k=%d, \\delta=%.2f', n_u, conf.k, conf.delta));
    
    legend('Location', 'bestoutside', 'NumColumns', 2);
    
    ylim([0 max(AEoL_NP_Sim(1,:))*1.3]); 
    
    drawnow;
end
