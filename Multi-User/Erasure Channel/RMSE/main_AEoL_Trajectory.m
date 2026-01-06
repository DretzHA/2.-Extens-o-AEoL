clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES
%  ========================================================================
conf.k = 16;          
conf.delta = 0.8;    
conf.velocities = 2 * ones(2, 1); % Exemplo: 2 m/s para 5 usuários

conf.Symbol_Duration = 1; 

% Parâmetros do Canal
avg_service_time = conf.k / (1 - conf.delta);
avg_service_time_sec = avg_service_time * conf.Symbol_Duration;
mu_eff = 1 / avg_service_time_sec; 

conf.rho_total_vec = 0.1:0.5:5; 
conf.lambda_total_vec = conf.rho_total_vec * mu_eff; 

conf.N_updates = 200;   
conf.N_MC = 3000; % Monte Carlo

% --- Configuração do Erro de Posicionamento (RMSE) ---
% Aqui definimos o desvio padrão do erro do sensor (em metros).
% O erro real de cada pacote será uma amostra de N(0, sigma^2).
conf.sensor_sigma = 3.0; 

%% ========================================================================
%  SIMULAÇÃO
%  ========================================================================
user_scenarios = [2]; 

for i = 1:length(user_scenarios)
    n_u = user_scenarios(i);
    conf.num_users = n_u;
    
    % Ajustar vetor de velocidades se necessário
    if length(conf.velocities) < n_u
        conf.velocities = 2 * ones(n_u, 1); 
    end
    
    fprintf('Simulating AEoL (Trajectory Based) for N=%d users...\n', n_u);
   
    % Chamada da função principal
    [AEoL_NP, AEoL_P, Sys_NP, Sys_P] = sim_mc_aeol_trajectory(conf);
    
    %% ====================================================================
    %  PLOTAGEM
    %  ====================================================================
    figure('Color', 'w'); hold on; grid on; box on;
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');
    
    % Cores
    c_np = [0.8500 0.3250 0.0980];
    c_p  = [0 0.4470 0.7410];
    
    % Plot Sistema
    plot(conf.lambda_total_vec, Sys_NP, '-o', 'Color', c_np, 'LineWidth', 2, 'DisplayName', 'System NP');
    plot(conf.lambda_total_vec, Sys_P,  '--s', 'Color', c_p,  'LineWidth', 2, 'DisplayName', 'System P');
    
    xlabel('Total Update Rate $\lambda_{tot}$', 'Interpreter', 'latex');
    ylabel('Average AEoL (m)', 'Interpreter', 'latex');
    title(sprintf('AEoL via Trajectory Error Integration (N=%d)', n_u), 'Interpreter', 'latex');
    legend('Location', 'best');
    
    drawnow;
end