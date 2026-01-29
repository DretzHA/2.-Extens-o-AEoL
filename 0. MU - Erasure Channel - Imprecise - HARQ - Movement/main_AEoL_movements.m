clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES DO CENÁRIO 2D (RWP)
%  ========================================================================
conf.k = 1;          % Quantidade de simbolos
conf.delta = 0.2;    % Probabilidade de erasure
conf.velocity = 1;    % Velocidade (m/s)
conf.area_side = 50; % Área 10x10m

conf.Symbol_Duration = 1; 

% Parâmetros do Canal (M/G/1)
avg_service_time = conf.k / (1 - conf.delta); 
avg_service_time_sec = avg_service_time * conf.Symbol_Duration; 
mu_eff = 1 / avg_service_time_sec; 
conf.mu = mu_eff;

conf.rho_total_vec = 0.5:0.5:10; 
conf.lambda_total_vec = conf.rho_total_vec * mu_eff; 

% --- Simulação ---
conf.Sim_Time = 50;   
conf.N_MC = 200;       % Monte Carlo
conf.sensor_sigma = 2; % Erro do sensor (m)
user_scenarios = [1]; % Cenários de teste
  
for i = 1:length(user_scenarios)
    n_u = user_scenarios(i); 
    conf.num_users = n_u;  
    
    fprintf('Simulando AEoL 2D (RWP) para N=%d usuários...\n', n_u);
   
    [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P, ...
     Sim_User_Counts_NP, Sim_User_Counts_P] = sim_mc_aeol_movement(conf);
    
    %% ====================================================================
    %  PLOTS
    %  ====================================================================
    figure('Name', sprintf('AEoL System 2D - %d Users', n_u), 'Color', 'w');
    subplot(2, 1, 1); hold on; grid on; box on;
    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    
    % Cores para cada usuário
    user_colors = lines(n_u);
    color_sys_p  = [0 0.4470 0.7410];      % Azul
    color_sys_np = [0.8500 0.3250 0.0980]; % Laranja
    
    % --- PLOT INDIVIDUAL (Background) ---
    for u = 1:n_u
       col_u = user_colors(u, :);
       
       % User NP (Linha fina, transparente)
       plot(conf.lambda_total_vec, Sim_AEoL_NP(u, :), ...
           'LineStyle', '-', 'Color', [col_u 0.3], 'LineWidth', 1.0, ... 
           'HandleVisibility', 'off'); 
           
       % User P (Linha pontilhada, transparente)
       plot(conf.lambda_total_vec, Sim_AEoL_P(u, :), ...
           'LineStyle', ':', 'Color', [col_u 0.6], 'LineWidth', 1.0, ...
           'HandleVisibility', 'off'); 
    end

    % --- PLOT SISTEMA (Destaque) ---
    % M/G/1 LCFS-NP System
    plot(conf.lambda_total_vec, Sim_Sys_NP, '-o', 'Color', color_sys_np, ...
        'MarkerFaceColor', color_sys_np, 'LineWidth', 2, 'DisplayName', 'Sys LCFS-NP');
    
    % M/G/1 LCFS-P System
    plot(conf.lambda_total_vec, Sim_Sys_P,  '--s', 'Color', color_sys_p, ...
        'MarkerFaceColor', color_sys_p, 'LineWidth', 2, 'DisplayName', 'Sys LCFS-P');
    
    xlabel('Taxa Total $\lambda_{tot}$', 'Interpreter', 'latex');
    ylabel('Age of Positioning (MSE Area)', 'Interpreter', 'latex');
    title(sprintf('AEoL: System vs Individual Users (N=%d)', n_u), 'Interpreter', 'latex');
    legend('Location', 'best');
    
    % Plot Throughput (Opcional, para validação)
    subplot(2, 1, 2); hold on; grid on; box on;
    set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex');
    plot(conf.lambda_total_vec, sum(Sim_User_Counts_NP,1)/conf.Sim_Time, '-o', 'DisplayName', 'Throughput NP');
    plot(conf.lambda_total_vec, sum(Sim_User_Counts_P,1)/conf.Sim_Time, '--s', 'DisplayName', 'Throughput P');
    xlabel('$\lambda_{tot}$', 'Interpreter', 'latex');
    ylabel('Packets/sec', 'Interpreter', 'latex');
    legend('Location', 'best');
    
    drawnow;
end