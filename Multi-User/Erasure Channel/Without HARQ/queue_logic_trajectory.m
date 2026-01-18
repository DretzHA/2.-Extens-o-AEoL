function [avg_aeol_user, avg_aeol_system, num_processed_total, actuation_time_per_user, processed_counts_per_user] = ...
    queue_logic_trajectory(arrivals, services, errors, velocities, type, T_max, is_successful)
    
    num_users = length(velocities);
    
    % Data cols: 1:Arr, 2:User, 3:Serv, 4:Err, 5:Success
    data_full = [arrivals, services, errors, is_successful];

    % Rodar Fila (O comportamento da fila não muda, a competição por recursos continua)
    if strcmp(type, 'NP')
        completed = run_lcfs_np(data_full, T_max);
    else
        completed = run_lcfs_p(data_full, T_max);
    end
    
    % completed columns: [FinishTime, GenTime, UserID, MeasError]
    num_processed_total = size(completed, 1);
    
    % Inicializar métricas
    actuation_time_per_user = zeros(num_users, 1);
    processed_counts_per_user = zeros(num_users, 1);
    avg_aeol_user = zeros(num_users, 1);
    
    % Velocidade do veículo (assumindo que todos rastreiam o mesmo alvo com mesma dinâmica)
    % Se sensores tivessem visões diferentes, seria mais complexo. 
    % Aqui usamos a velocidade do alvo (v_mean).
    v_target = mean(velocities); 

    % ---------------------------------------------------------------------
    % 1. Calcular AEoL Individual (Como se cada sensor estivesse sozinho)
    % ---------------------------------------------------------------------
    for u = 1:num_users
        % --- TRATAMENTO DE ERRO 2: Completed Vazio ---
        % Se completed for [], tentar acessar completed(:,3) quebra o código.
        if isempty(completed)
            pkts_u = [];
        else
            % Só faz a indexação se tiver dados
            idx_u = (completed(:, 3) == u);
            pkts_u = completed(idx_u, :);
        end
        
        processed_counts_per_user(u) = size(pkts_u, 1);
        
        if isempty(pkts_u)
            % Se o sensor não entregou nada, erro máximo
            avg_aeol_user(u) = 0.5 * v_target * T_max;
        else
            % Ordena e calcula área
            pkts_u = sortrows(pkts_u, 1);
            area_u = calc_position_error_area(pkts_u, T_max, v_target);
            avg_aeol_user(u) = area_u / T_max;
        end
    end
    
    % ---------------------------------------------------------------------
    % 2. Calcular AEoL do SISTEMA (Fusão de Sensores)
    % ---------------------------------------------------------------------
    % AQUI ESTÁ A CORREÇÃO: Usamos 'completed' inteiro (todos os sensores)
    % para traçar a trajetória real do veículo único.
    
    if isempty(completed)
        avg_aeol_system = 0.5 * v_target * T_max;
    else
        % Ordena todas as atualizações cronologicamente (não importa a origem)
        all_updates_sorted = sortrows(completed, 1);
        
        % Calcula a área combinada
        area_system = calc_position_error_area(all_updates_sorted, T_max, v_target);
        avg_aeol_system = area_system / T_max;
    end
end

% --- Funções Auxiliares (calc_position_error_area, run_lcfs_np, etc) ---
% MANTENHA AS FUNÇÕES ABAIXO IGUAIS ÀS QUE ENVIEI ANTERIORMENTE
% (Elas já suportam o cálculo correto de área se passarmos a lista certa)

function area = calc_position_error_area(updates, T_max, v)
    % updates: [FinTime, GenTime, UserID, MeasError]
    t_prev = 0;
    last_gen = 0; 
    last_err = 0; 
    area = 0;
    
    n_up = size(updates, 1);
    
    for k = 1:n_up
        t_curr = updates(k, 1); 
        t_gen_new = updates(k, 2);
        err_new   = updates(k, 4);
        
        slope = v;
        offset = -v * last_gen + last_err;
        
        seg_area = integrate_abs_linear(slope, offset, t_prev, t_curr);
        area = area + seg_area;
        
        t_prev = t_curr;
        last_gen = t_gen_new;
        last_err = err_new;
    end
    
    slope = v;
    offset = -v * last_gen + last_err;
    area = area + integrate_abs_linear(slope, offset, t_prev, T_max);
end

function val = integrate_abs_linear(a, b, t1, t2)
    t_z = -b/a;
    if t_z > t1 && t_z < t2
        val = integral_simple(a, b, t1, t_z) + integral_simple(a, b, t_z, t2);
    else
        val = integral_simple(a, b, t1, t2);
    end
end

function val = integral_simple(a, b, t_s, t_e)
    t_m = (t_s + t_e)/2; 
    val_m = a*t_m + b;
    if val_m < 0
        F = @(t) -1 * (0.5*a*t^2 + b*t);
    else
        F = @(t) (0.5*a*t^2 + b*t);
    end
    val = F(t_e) - F(t_s);
end

% --- Lógicas de Fila (Reutilizar as versões com IsSuccess da resposta anterior) ---
function completed = run_lcfs_np(data, T_max)
    N = size(data, 1); completed = []; served_mask = false(N, 1); time_now = 0;
    while true
        if time_now >= T_max, break; end
        queue = find(data(:,1) <= time_now & ~served_mask);
        if isempty(queue)
            upcoming = find(data(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming), break; end
            time_now = data(upcoming, 1); continue;
        end
        target = queue(end);
        if length(queue) > 1, served_mask(queue(1:end-1)) = true; end
        fin_t = time_now + data(target, 3);
        if fin_t <= T_max
            if data(target, 5) == 1 % Check Success
                completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)];
            end
            time_now = fin_t;
        else
            time_now = fin_t; 
        end
        served_mask(target) = true;
    end
end

function completed = run_lcfs_p(data, T_max)
    N = size(data, 1); completed = []; 
    if N==0, return; end
    for k = 1:N-1
        ft = data(k,1) + data(k,3);
        if ft <= data(k+1,1) && ft <= T_max
            if data(k, 5) == 1
                completed = [completed; ft, data(k,1), data(k,2), data(k,4)];
            end
        end
    end
    ft = data(N,1) + data(N,3);
    if ft <= T_max
        if data(N, 5) == 1
            completed = [completed; ft, data(N,1), data(N,2), data(N,4)];
        end
    end
end