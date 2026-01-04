%% main.m
clc;
clear;
close all;

%% ========================================================================
%  1. CONFIGURAÇÕES (PARÂMETROS)
%  ========================================================================
mu = 1;                     % Taxa de serviço (pacotes/s)
velocities = [1];           % Velocidade do alvo (m/s)
rho_vec = 0.05:0.05:1.2;    % Taxa de utilização (rho = lambda/mu)
lambda_vec = rho_vec * mu;  % Taxa de chegada

N_updates = 200;           % Número de pacotes por rodada (para estabilidade)
N_MC = 500;                 % Número de iterações Monte Carlo

%% ========================================================================
%  2. CÁLCULOS
%  ========================================================================

% 2.1. Executar Cálculos Analíticos (Teoria)
fprintf('Calculando resultados teóricos...\n');
[AEoL_woP_Analytic, AEoL_withP_Analytic] = ...
    compute_analytical_results(velocities, mu, rho_vec, lambda_vec);

% 2.2. Executar Simulação (Monte Carlo)
fprintf('Executando Monte Carlo (%d iterações)...\n', N_MC);
[AEoL_woP_Sim, AEoL_withP_Sim] = ...
    run_monte_carlo_simulation(velocities, mu, lambda_vec, N_MC, N_updates);

fprintf('Processamento concluído.\n');

%% ========================================================================
%  3. PLOTAGEM
%  ========================================================================
custom_colors = {'#D95319', '#7E2F8E'}; % Laranja, Roxo

figure('Name', 'AEoL Modular', 'Color', 'w');
hold on;
set(gca, 'FontSize', 12);

for i = 1:length(velocities)
    v = velocities(i);
    
    % Sem Preempção
    plot(lambda_vec, AEoL_woP_Analytic(i, :), '-', 'LineWidth', 1.5, ...
        'Color', custom_colors{1}, 'DisplayName', 'Analyt. (Sem Preempção)');
    plot(lambda_vec, AEoL_woP_Sim(i, :), '*', 'LineWidth', 1.5, ...
        'Color', custom_colors{1}, 'DisplayName', 'Sim. (Sem Preempção)');
        
    % Com Preempção
    plot(lambda_vec, AEoL_withP_Analytic(i, :), '--', 'LineWidth', 1.5, ...
        'Color', custom_colors{2}, 'DisplayName', 'Analyt. (Com Preempção)');
    plot(lambda_vec, AEoL_withP_Sim(i, :), 'o', 'LineWidth', 1.5, ...
        'Color', custom_colors{2}, 'DisplayName', 'Sim. (Com Preempção)');
end

xlabel('Taxa de Chegada \lambda');
ylabel('AEoL [m]');
title('AEoL: Comparação Analítico vs Simulação');
legend('show', 'Location', 'north');
grid on; box on;
hold off;