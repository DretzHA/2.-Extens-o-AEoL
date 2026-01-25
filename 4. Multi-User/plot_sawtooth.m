%% plot_sawtooth_multi_visual.m
% Visualização Dente de Serra (AoI) para 2 Usuários Simultâneos
% Compara LCFS-NP (Sem Preempção) vs LCFS-P (Com Preempção)
clc; clear; close all;

%% 1. CONFIGURAÇÕES
conf.mu = 1.0;             % Taxa de serviço
conf.num_users = 3;        % 2 Usuários para visualização clara
conf.rho = 0.5;           % Carga alta para forçar interação/disputa
conf.N_updates = 15;       % Poucos pacotes para dar "zoom" no início
conf.t_limit = 12;         % Limite de tempo no eixo X para o plot

% Cores para os usuários (Azul e Laranja padrão MATLAB)
conf.colors = [0, 0.4470, 0.7410;  % User 1 (Azul)
               0.8500, 0.3250, 0.0980]; % User 2 (Laranja)

rng(42); % Semente fixa para garantir que o gráfico seja interessante

%% 2. GERAÇÃO DE TRÁFEGO
lambda_tot = conf.rho * conf.mu;
lambda_user = lambda_tot / conf.num_users;

all_arrivals = [];
for u = 1:conf.num_users
    dt = exprnd(1/lambda_user, conf.N_updates, 1);
    t_arr = cumsum(dt);
    all_arrivals = [all_arrivals; t_arr, ones(conf.N_updates, 1)*u];
end

% Ordena cronologicamente (Fila Única Compartilhada)
all_arrivals = sortrows(all_arrivals, 1);

% Gera serviços
n_pkts = size(all_arrivals, 1);
services = exprnd(1/conf.mu, n_pkts, 1);

%% 3. SIMULAÇÃO (Lógica de Filas)

% NP (Non-Preemptive)
completed_NP = run_lcfs_np(all_arrivals, services);

% P (Preemptive)
completed_P = run_lcfs_p(all_arrivals, services);

%% 4. PLOTAGEM GERAL
figure('Name', 'Multi-User AoI Dynamics', 'Color', 'w', 'Position', [100 50 1000 800]);

% --- SUBPLOT 1: SEM PREEMPÇÃO ---
subplot(2,1,1);
hold on;
title('LCFS Non-Preemptive');
% Loop para desenhar cada usuário
for u = 1:conf.num_users
    draw_user_layer(all_arrivals, completed_NP, u, conf.colors(u,:), conf.t_limit);
end
setup_axes(conf.t_limit);

% --- SUBPLOT 2: COM PREEMPÇÃO ---
subplot(2,1,2);
hold on;
title('LCFS Preemptive');
% Loop para desenhar cada usuário
for u = 1:conf.num_users
    draw_user_layer(all_arrivals, completed_P, u, conf.colors(u,:), conf.t_limit);
end
setup_axes(conf.t_limit);

fprintf('Gráfico gerado com sucesso.\n');

%% ========================================================================
%  FUNÇÕES AUXILIARES DE PLOTAGEM
%  =======================================================================
function draw_user_layer(arrivals, completed, u_id, color_code, t_lim)
    % Filtra chegadas deste usuário
    my_arr = arrivals(arrivals(:,2) == u_id, 1);
    
    % Filtra pacotes completados deste usuário: [Finish, Gen, UserID]
    mask = (completed(:,3) == u_id);
    my_comp = completed(mask, 1:2); 
    my_comp = sortrows(my_comp, 1); % Ordena por tempo de saída
    
    % --- 1. Plotar Chegadas (Triângulos no eixo X) ---
    % Pequeno offset vertical para não sobrepor se chegarem juntos
    offset = (u_id - 1) * -0.2; 
    stem(my_arr, ones(size(my_arr))*0.2 + offset, 'Marker', '^', ...
        'Color', color_code, 'LineStyle', 'none', 'MarkerFaceColor', 'w', ...
        'HandleVisibility', 'off'); % Oculta da legenda para não poluir
    
    % --- 2. Construir Curva AoI ---
    X = [0];
    Y = [0];
    t_prev = 0;
    last_gen = 0;
    
    for i = 1:size(my_comp, 1)
        t_finish = my_comp(i, 1);
        t_gen    = my_comp(i, 2);
        
        % Só consideramos pacotes "frescos" (geração posterior à última conhecida)
        if t_gen > last_gen
            % Sobe até o momento da entrega
            X = [X, t_finish];
            Y = [Y, (t_finish - last_gen)];
            
            % Cai para a nova idade
            new_aoi = t_finish - t_gen;
            X = [X, t_finish];
            Y = [Y, new_aoi];
            
            % Marcador de entrega (Bolinha cheia)
            plot(t_finish, new_aoi, 'o', 'MarkerFaceColor', color_code, ...
                'MarkerEdgeColor', 'none', 'MarkerSize', 6, 'HandleVisibility', 'off');
            
            last_gen = t_gen;
            t_prev = t_finish;
        end
    end
    
    % Linha final até o limite
    X = [X, t_lim];
    Y = [Y, t_lim - last_gen];
    
    % Plota a linha grossa
    plot(X, Y, '-', 'Color', color_code, 'LineWidth', 2, ...
        'DisplayName', sprintf('User %d', u_id));
end

function setup_axes(t_lim)
    grid on;
    xlabel('T (s)');
    ylabel('AoI');
    xlim([0 t_lim]);
    % Limite Y dinâmico mas controlado
    ylim([0 t_lim/1.5]); 
    legend('show', 'Location', 'northwest');
end

%% ========================================================================
%  LÓGICA DE FILAS (Core Simulation)
%  ========================================================================
function completed = run_lcfs_np(arr, srv)
    N = size(arr, 1);
    completed = []; 
    time_now = 0;
    server_free = 0;
    served_mask = false(N, 1);
    processed_count = 0;
    
    while processed_count < N
        if time_now < server_free
            time_now = server_free;
        end
        queue_idx = find(arr(:,1) <= time_now & ~served_mask);
        if isempty(queue_idx)
            next_idx = find(~served_mask, 1, 'first');
            if isempty(next_idx), break; end
            time_now = arr(next_idx, 1);
        else
            target = queue_idx(end); % LCFS
            start_t = time_now;
            finish_t = start_t + srv(target);
            completed = [completed; finish_t, arr(target, 1), arr(target, 2)];
            served_mask(target) = true;
            server_free = finish_t;
            time_now = finish_t;
            processed_count = processed_count + 1;
        end
    end
end

function completed = run_lcfs_p(arr, srv)
    N = size(arr, 1);
    completed = [];
    % Preempção global: Só completa se finish < próxima chegada QUALQUER
    for k = 1:(N-1)
        arr_curr = arr(k, 1);
        arr_next = arr(k+1, 1);
        srv_time = srv(k);
        pot_fin = arr_curr + srv_time;
        
        if pot_fin <= arr_next
            completed = [completed; pot_fin, arr_curr, arr(k, 2)];
        end
    end
    % Último pacote
    completed = [completed; arr(N,1)+srv(N), arr(N,1), arr(N,2)];
end