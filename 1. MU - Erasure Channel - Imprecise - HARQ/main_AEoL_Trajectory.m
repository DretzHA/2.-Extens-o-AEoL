clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES
%  ========================================================================
conf.k = 16;      %quantidade de simbolos na informação    
conf.delta = 0.20;    %probabilidade de erasure
velocity = 1; %velocidade dos usuários
 
conf.Symbol_Duration = 1e-1; %duração de um simbolo

% Parâmetros do Canal
avg_service_time = conf.k / (1 - conf.delta); %tempo médio de serviço (Eq. 4 - Yates 2017)
avg_service_time_sec = avg_service_time * conf.Symbol_Duration; %tempo em s
mu_eff = 1 / avg_service_time_sec; %Taxa de serviço
conf.mu = mu_eff;

conf.rho_total_vec = 0.1:0.25:5; % carga do servidor
conf.lambda_total_vec = conf.rho_total_vec * mu_eff; %taxa de atualizações

% --- Simulação ---
conf.Sim_Time = 500;   %Tempo de simulação
conf.N_MC = 500;      %Número de Simulações Monte Carlo
conf.sensor_sigma = 0; %Desvio Padrão do Erro de Localização
user_scenarios = [2]; %Total de sources

for i = 1:length(user_scenarios)
    n_u = user_scenarios(i); %número de usuários atuais
    conf.num_users = n_u;  
    conf.velocities = velocity * ones(n_u, 1); 

    %plot_sawtooth_trajectory(conf, n_u);
    %plot_trajectory_tracking(conf);
    
    fprintf('Simulating AEoL for N=%d users\n', n_u, conf.Sim_Time);
   
    % Chama função para simulação
    [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P, Thr_NP, Thr_P, ...
     Sim_User_Counts_NP, Sim_User_Counts_P] = sim_mc_aeol_trajectory(conf);
    
    % Calcular totais absolutos
    Total_Pkts_NP = Thr_NP * conf.Sim_Time;
    Total_Pkts_P  = Thr_P  * conf.Sim_Time;

    %% ---------------------------------------------------------------------
    %  THEORY CALCULATION
    % ---------------------------------------------------------------------
    aeol_theory_np = zeros(size(conf.lambda_total_vec));
    aeol_theory_p  = zeros(size(conf.lambda_total_vec));
    v = conf.velocities(1);
    
    for j = 1:length(conf.lambda_total_vec)
        lambda_val = conf.lambda_total_vec(j); 

        %%Teste Teorema 2 - Yates 2019
        % lambda_val_user = conf.lambda_total_vec(j)/n_u; 
        % 
        % rho_val = conf.rho_total_vec(j);
        % rho_val_user = conf.rho_total_vec(j)/n_u; 
        
        lambda_per_slot = lambda_val * conf.Symbol_Duration;
        
        % --- NP Theory (Slots) --- Eq 5 - Najm 2019
        num = lambda_per_slot * conf.k * (conf.k + conf.delta);
        denom = 2 * (1 - conf.delta) * (lambda_per_slot * conf.k + 1 - conf.delta);
        val_np_slots = (1/lambda_per_slot) + (conf.k/(1-conf.delta)) + (num/denom);
        aeol_theory_np(j) = v * val_np_slots * conf.Symbol_Duration;
        
        % --- P Theory (Slots) --- Eq 17 - Najm 2019
        term_base = (exp(lambda_per_slot) - conf.delta) / (1 - conf.delta);
        val_p_slots = (1/lambda_per_slot) * (term_base ^ conf.k);
        aeol_theory_p(j) = v * val_p_slots * conf.Symbol_Duration;
    
        %%Teste Teorema 2 - Yates 2019
        % val_p_slots = (1/mu_eff) * (1+rho_val)*(1/rho_val_user);
        % aeol_theory_p(j) = v * val_p_slots * conf.Symbol_Duration;


    end

    %% ====================================================================
    %  FIGURAS
    %  ====================================================================
    fig_name_aeol = sprintf('AEoL Results - %d Users', n_u);
    figure('Name', fig_name_aeol, 'Color', 'w', 'Position', [50 100 800 600]);
    
    hold on; grid on; grid minor; box on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.1, 'TickLabelInterpreter', 'latex');
    
    color_sys_p  = [0 0.4470 0.7410];      % Azul
    color_sys_np = [0.8500 0.3250 0.0980]; % Laranja
    color_theo   = [0.15 0.15 0.15];       % Cinza escuro
    
    user_colors = lines(n_u); 
    
    % --- Theory Plots ---
    plot(conf.lambda_total_vec, aeol_theory_np, ...
        'LineStyle', '-', 'Color', color_theo, 'LineWidth', 2.5, ...
        'DisplayName', 'Theory (NP)');
        
    plot(conf.lambda_total_vec, aeol_theory_p, ...
        'LineStyle', '--', 'Color', color_theo, 'LineWidth', 2.5, ...
        'DisplayName', 'Theory (P)');
    
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
    hold on; grid on; box on;
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