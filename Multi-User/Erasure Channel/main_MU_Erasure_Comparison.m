clc; clear; close all;

%% ========================================================================
%  CONFIGURAÇÕES GERAIS
%  ========================================================================
conf.k = 16;          
conf.delta = 0.10;    
conf.velocities = 1*ones(10, 1); 

conf.Symbol_Duration = 0.10; % Ex: 1 ms para transmitir 1 símbolo (ou pacote base)

% mu = service rate | rho = system load | lambda = arrival rate
% Eq. 4 - Timely Updates over an Erasure Channel
avg_service_time = conf.k / (1 - conf.delta);
avg_service_time_sec = avg_service_time * conf.Symbol_Duration;

mu_eff = 1 / avg_service_time_sec; 

conf.rho_total_vec = 0.1:0.5:5; 
conf.lambda_total_vec = conf.rho_total_vec * mu_eff; 

conf.N_updates = 100;   
conf.N_MC = 10000; % Reduzi levemente para teste rápido, aumente se necessário

user_scenarios = [1 2 3]; 

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
   
    % Roda a simulação Monte Carlo (necessita funções atualizadas anteriormente)
    [AEoL_NP_Sim, AEoL_P_Sim, AEoL_Sys_NP, AEoL_Sys_P] = sim_mc_multi_erasure(conf);
    
    % ---------------------------------------------------------------------
    % CÁLCULO TEÓRICO (Escalado pelo tempo de símbolo)
    % ---------------------------------------------------------------------
    aeol_theory_np = zeros(size(conf.lambda_total_vec));
    aeol_theory_p  = zeros(size(conf.lambda_total_vec));
    
    for j = 1:length(conf.lambda_total_vec)
        lambda_val = conf.lambda_total_vec(j); 
        
        % Converte lambda para (1/slot) para usar fórmulas discretas
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
    %  GERAÇÃO DA FIGURA COMBINADA (ESTILO PROFISSIONAL)
    %  ===================================================================
    fig_name = sprintf('Combined Results (NP vs P) - %d Users', n_u);
    % Aumentei um pouco a largura para a legenda caber melhor na direita
    figure('Name', fig_name, 'Color', 'w', 'Position', [100 100 1200 700]);
    
    % Configuração de eixos mais limpa
    hold on; 
    grid on; grid minor; box on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.1, 'TickLabelInterpreter', 'latex');
    
    % --- DEFINIÇÃO DE CORES (Paleta "Publication Safe") ---
    % Azul Escuro para P (System) | Vermelho Tijolo para NP (System)
    color_sys_p  = [0 0.4470 0.7410];      
    color_sys_np = [0.8500 0.3250 0.0980]; 
    color_theo   = [0.15 0.15 0.15];       % Cinza Quase Preto
    
    % Gera cores para os usuários (baseadas em 'lines' mas com transparência simulada se fossem muitos)
    user_colors = lines(n_u); 
    
    % --- 1. Plotar Curvas Teóricas (Camada de Fundo) ---
    % Theory NP: Linha Sólida Cinza Escuro
    plot(conf.lambda_total_vec, aeol_theory_np, ...
        'LineStyle', '-', 'Color', color_theo, 'LineWidth', 2.5, ...
        'DisplayName', 'Theory (NP)');
        
    % Theory P: Linha Tracejada Cinza Escuro
    plot(conf.lambda_total_vec, aeol_theory_p, ...
        'LineStyle', '--', 'Color', color_theo, 'LineWidth', 2.5, ...
        'DisplayName', 'Theory (P)');

    % --- 2. Plotar Curvas dos Usuários (Camada do Meio) ---
    % Linhas finas para não ofuscar o sistema
    for u = 1:n_u
        col_u = user_colors(u, :);
        % User NP (Linha fina contínua)
        plot(conf.lambda_total_vec, AEoL_NP_Sim(u, :), ...
            'LineStyle', '-', 'Color', [col_u 0.3], 'LineWidth', 0.8, ... % 0.3 é transparência (alpha) em versões novas
            'HandleVisibility', 'off'); 
        % User P (Linha fina pontilhada)
        plot(conf.lambda_total_vec, AEoL_P_Sim(u, :), ...
            'LineStyle', ':', 'Color', [col_u 0.6], 'LineWidth', 1.0, ...
            'HandleVisibility', 'off'); 
    end
    
    % "Dummy plots" para Legenda dos Usuários (em cinza para representar a classe)
    plot(NaN, NaN, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'DisplayName', 'Indiv. Users (NP)');
    plot(NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'DisplayName', 'Indiv. Users (P)');

    % --- 3. Plotar Curvas do SISTEMA (Destaque Principal) ---
    % System NP: Vermelho, Solido, Marcador Circular Cheio
    plot(conf.lambda_total_vec, AEoL_Sys_NP, ...
        'LineStyle', '-', 'Color', color_sys_np, 'LineWidth', 2.2, ...
        'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', color_sys_np, ...
        'DisplayName', 'SYSTEM (NP)');

    % System P: Azul, Tracejado, Marcador Quadrado Cheio
    plot(conf.lambda_total_vec, AEoL_Sys_P, ...
        'LineStyle', '--', 'Color', color_sys_p, 'LineWidth', 2.2, ...
        'Marker', 's', 'MarkerSize', 6, 'MarkerFaceColor', color_sys_p, ...
        'DisplayName', 'SYSTEM (P)');

    % --- Formatação Final ---
    xlabel('Total Update Rate $\lambda_{tot}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Average AEoL $\bar{\Delta}$', 'Interpreter', 'latex', 'FontSize', 14);
    
    title_str = sprintf('Scenario: %d Users (Shared Object) | k=%d, $\\delta=%.1f$', ...
                  n_u, conf.k, conf.delta);
    title(title_str, 'Interpreter', 'latex', 'FontSize', 14);
    
    % Ajuste da Legenda
    % 'NumColumns', 2 ajuda a economizar altura
    lgd = legend('Location', 'northeast', 'NumColumns', 2, 'FontSize', 10);
    lgd.ItemTokenSize = [20, 18]; % Reduz o tamanho da linha na legenda para ficar mais compacto
    
    % Limites
    yl = ylim;
    ylim([0 yl(2)*1.05]); % Dá um respiro de 5% no topo
    xlim([min(conf.lambda_total_vec) max(conf.lambda_total_vec)]);
    
    drawnow;
end