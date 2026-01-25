clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES
%  ========================================================================
conf.k = 16; %quantidade de símbolos na informação         
conf.delta = 0.0;     %probabilidade de erasure
velocity = 1;  %velocidade dos usuários

conf.Symbol_Duration = 1e-1; %duranção de cada símbolo

% Parâmetros do Canal
avg_service_time_sec = conf.k * conf.Symbol_Duration;
mu_eff = 1 / avg_service_time_sec; 

conf.rho_total_vec = 0.1:0.25:5; 
conf.lambda_total_vec = conf.rho_total_vec * mu_eff; 

% --- Simulação baseada em TEMPO ---
conf.Sim_Time = 500;   
conf.N_MC = 500;      
conf.sensor_sigma = 0; 



%% ========================================================================
%  SIMULAÇÃO
%  ========================================================================
user_scenarios = [2]; 

for i = 1:length(user_scenarios)
    n_u = user_scenarios(i);
    conf.num_users = n_u;  
    conf.velocities = velocity * ones(n_u, 1); 

    % plot_sawtooth_woHARQ(conf, n_u);
    % plot_trajectory_woHARQ(conf);
    
    fprintf('Simulating AEoL for N=%d users, Time=%.1fs...\n', n_u, conf.Sim_Time);
   
    % Chama simulação (com todas as saídas definidas anteriormente)
    [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P, Thr_NP, Thr_P, ...
     Sim_Act_Time_NP, Sim_Act_Time_P, Sim_User_Counts_NP, Sim_User_Counts_P] = sim_mc_aeol_trajectorywo_HARQ(conf);
    
    % Calcular totais absolutos
    Total_Pkts_NP = Thr_NP * conf.Sim_Time;
    Total_Pkts_P  = Thr_P  * conf.Sim_Time;

    %% ====================================================================
    %  FIGURA 1: Resultados de AEoL (Sistema + Usuários + Teoria)
    %  ====================================================================
    fig_name_aeol = sprintf('AEoL Results - %d Users', n_u);
    figure('Name', fig_name_aeol, 'Color', 'w', 'Position', [50 100 800 600]);
    
    hold on; grid on; grid minor; box on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.1, 'TickLabelInterpreter', 'latex');
    
    color_sys_p  = [0 0.4470 0.7410];      % Azul
    color_sys_np = [0.8500 0.3250 0.0980]; % Laranja
    
    user_colors = lines(n_u);   
    
    % --- Plot: Usuários Individuais ---
    for u = 1:n_u
        col_u = user_colors(u, :);
        
        % User NP (Linha sólida, transparente)
        plot(conf.lambda_total_vec, Sim_AEoL_NP(u, :), ...
            'LineStyle', '-', 'Color', [col_u 0.3], 'LineWidth', 1.0, ... 
            'HandleVisibility', 'off'); 
            
        % User P (Linha pontilhada, transparente)
        plot(conf.lambda_total_vec, Sim_AEoL_P(u, :), ...
            'LineStyle', ':', 'Color', [col_u 0.6], 'LineWidth', 1.0, ...
            'HandleVisibility', 'off'); 
    end
    
    % "Dummy plots" para criar a legenda cinza dos usuários
    plot(NaN, NaN, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'DisplayName', 'Indiv. Users (NP)');
    plot(NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'DisplayName', 'Indiv. Users (P)');
    
    % --- Plot: Sistema (Média) ---
    % System NP
    plot(conf.lambda_total_vec, Sim_Sys_NP, ...
        'LineStyle', '-', 'Color', color_sys_np, 'LineWidth', 2.5, ...
        'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', color_sys_np, ...
        'DisplayName', 'SYSTEM (NP)');
        
    % System P
    plot(conf.lambda_total_vec, Sim_Sys_P, ...
        'LineStyle', '--', 'Color', color_sys_p, 'LineWidth', 2.5, ...
        'Marker', 's', 'MarkerSize', 6, 'MarkerFaceColor', color_sys_p, ...
        'DisplayName', 'SYSTEM (P)');
    
    % Formatação
    xlabel('Total Update Rate', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Average AEoL (m)', 'Interpreter', 'latex', 'FontSize', 14);
    title(sprintf('AEoL Analysis: %d Users', n_u), 'Interpreter', 'latex', 'FontSize', 14);
    
    lgd = legend('Location', 'northeast', 'NumColumns', 1, 'FontSize', 10);
    
    % Ajustar limites
    yl = ylim;
    ylim([0 yl(2)*1.05]); 
    xlim([min(conf.lambda_total_vec) max(conf.lambda_total_vec)]);

    %% ====================================================================
    %  FIGURA 2: Throughput e Quantidade de Pacotes
    %  ====================================================================
    fig_name_thr = sprintf('Throughput Stats - %d Users', n_u);
    figure('Name', fig_name_thr, 'Color', 'w', 'Position', [900 100 600 600]);
    
    % --- Subplot 1: Vazão (Updates / sec) ---
    subplot(2, 1, 1); hold on; grid on; box on;
    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    
    plot(conf.lambda_total_vec, Thr_NP, '-o', 'Color', color_sys_np, 'LineWidth', 2, 'DisplayName', 'NP');
    plot(conf.lambda_total_vec, Thr_P,  '--s', 'Color', color_sys_p,  'LineWidth', 2, 'DisplayName', 'P');
    
    ylabel('Throughput (updates/sec)', 'Interpreter', 'latex');
    title('Effective Update Rate', 'Interpreter', 'latex');
    legend('Location', 'best');
    
    % --- Subplot 2: Total Absoluto de Pacotes Processados ---
    subplot(2, 1, 2); hold on; grid on; box on;
    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    
    plot(conf.lambda_total_vec, Total_Pkts_NP, '-o', 'Color', color_sys_np, 'LineWidth', 2, 'DisplayName', 'NP');
    plot(conf.lambda_total_vec, Total_Pkts_P,  '--s', 'Color', color_sys_p,  'LineWidth', 2, 'DisplayName', 'P');
    
    xlabel('Total Update Rate $\lambda_{tot}$', 'Interpreter', 'latex');
    ylabel('Total Updates Processed (Count)', 'Interpreter', 'latex');
    title(sprintf('Total Packets Processed (T = %.0fs)', conf.Sim_Time), 'Interpreter', 'latex');
    legend('Location', 'best');
    
    drawnow;

    %% ====================================================================
    %  FIGURA 3: Estatísticas por Usuário (Extra)
    %  ====================================================================
    figure('Name', sprintf('User Stats - %d Users', n_u), 'Color', 'w', 'Position', [150 50 1000 500]);
    
    subplot(1,2,1); hold on; grid on; box on;
    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    for u = 1:n_u
        col = user_colors(u, :);
        plot(conf.lambda_total_vec, Sim_Act_Time_NP(u, :), '-', 'Color', col, 'LineWidth', 1.5);
        plot(conf.lambda_total_vec, Sim_Act_Time_P(u, :), '--', 'Color', col, 'LineWidth', 1.5);
    end
    xlabel('$\lambda_{tot}$', 'Interpreter','latex');
    ylabel('Effective Actuation Time (s)', 'Interpreter','latex');
    title(sprintf('Actuation Duration (Total T=%.0fs)', conf.Sim_Time), 'Interpreter','latex');

    subplot(1,2,2); hold on; grid on; box on;
    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    for u = 1:n_u
        col = user_colors(u, :);
        plot(conf.lambda_total_vec, Sim_User_Counts_NP(u, :), '-', 'Color', col, 'LineWidth', 1.5);
        plot(conf.lambda_total_vec, Sim_User_Counts_P(u, :), '--', 'Color', col, 'LineWidth', 1.5);
    end
    xlabel('$\lambda_{tot}$', 'Interpreter','latex');
    ylabel('Processed Packets per User', 'Interpreter','latex');
    title('Packets per User', 'Interpreter','latex');
    
    drawnow;
end
