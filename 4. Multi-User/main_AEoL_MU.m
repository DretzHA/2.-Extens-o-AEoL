clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES
%  ========================================================================
conf.mu = 1;                       
conf.rho_total_vec = 0.1:0.1:1;
conf.lambda_total_vec = conf.rho_total_vec * conf.mu;

conf.N_updates = 100;              
conf.N_MC = 1; % Aumentei um pouco para estabilizar a simulação

user_scenarios = [2];
colors = lines(length(user_scenarios)); % Cores: Azul, Laranja, Amarelo

% Preparação da Figura
figure('Name', 'Comparação Completa AEoL', 'Color', 'w', 'Position', [100 100 1000 600]);
hold on;
set(gca, 'FontSize', 11);
grid on; box on;

%% ========================================================================
%  2. LOOP DE SIMULAÇÃO E PLOTAGEM
%  ========================================================================
fprintf('--- Iniciando Simulação Comparativa ---\n');

for i = 1:length(user_scenarios)
    n_u = user_scenarios(i);
    cor = colors(i, :);
    
    % Configura cenário atual
    conf.num_users = n_u;
    conf.velocities = ones(1, conf.num_users) * 1.0; 
    
    fprintf('Processando N=%d usuários...\n', n_u);
    
    % --- 2.1. Cálculos ---
    % Analítico
    [AEoL_NP_Analytic, AEoL_P_Analytic] = analytic_multi(conf);
    
    % Numérico (Simulação)
    % Certifique-se que a função 'sim_mc_multi' existe e retorna os vetores corretos
    [AEoL_NP_Sim, AEoL_P_Sim] = sim_mc_multi(conf);
    
    % Índice do usuário de referência (User 1)
    u_idx = 1; 
    
    % --- 2.2. Plotagem com Legendas Automáticas ---
    
    % GRUPO 1: SEM PREEMPÇÃO (Non-Preemptive)
    % Linha Analítica (Sólida)
    plot(conf.lambda_total_vec, AEoL_NP_Analytic(u_idx, :), ...
        'LineStyle', '-', 'Color', cor, 'LineWidth', 2, ...
        'DisplayName', sprintf('Users=%d NP (Analyt.)', n_u));
    
    % Pontos Simulados (x) - Sem linha conectando
    plot(conf.lambda_total_vec, AEoL_NP_Sim(u_idx, :), ...
        'Marker', 'x', 'LineStyle', 'none', 'Color', cor, 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'DisplayName', sprintf('Users=%d NP (Sim.)', n_u));
    
    % GRUPO 2: COM PREEMPÇÃO (Preemptive)
    % Linha Analítica (Tracejada)
    plot(conf.lambda_total_vec, AEoL_P_Analytic(u_idx, :), ...
        'LineStyle', '--', 'Color', cor, 'LineWidth', 2, ...
        'DisplayName', sprintf('Users=%d P (Analyt.)', n_u));
        
    % Pontos Simulados (Círculo) - Sem linha conectando
    plot(conf.lambda_total_vec, AEoL_P_Sim(u_idx, :), ...
        'Marker', 'o', 'LineStyle', 'none', 'Color', cor, 'LineWidth', 1.5, ...
        'MarkerSize', 6, 'DisplayName', sprintf('Users=%d P (Sim.)', n_u));
end

%% ========================================================================
%  3. FORMATAÇÃO FINAL
%  ========================================================================
xlabel('Update Rate (\lambda_{tot})');
ylabel('AEoL [m]');
title('AEoL LCFS: Analytical vs Simulated');

% A mágica da legenda organizada:
% 'NumColumns' distribui os itens lado a lado para não ficar uma lista vertical gigante
lgd = legend('show');
set(lgd, 'Location', 'northwest', 'NumColumns', 3, 'FontSize', 9);

ylim([0 max(AEoL_NP_Analytic(1,:))*1.2]); % Ajuste de zoom vertical se necessário

hold off;
fprintf('--- Fim ---\n');