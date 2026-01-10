function [avg_aeol_user, avg_aeol_system, num_processed_total, actuation_time_per_user, processed_counts_per_user] = queue_logic_trajectory(arrivals, services, errors, velocities, type, T_max)
    % Retorna:
    % - actuation_time_per_user: vetor [N_users x 1] com o tempo total que cada usuario ficou "ativo"
    % - processed_counts_per_user: vetor [N_users x 1] com contagem de pacotes
    
    num_users = length(velocities);
    data_full = [arrivals, services, errors];

    % Rodar Fila
    if strcmp(type, 'NP')
        completed = run_lcfs_np(data_full, T_max);
    else
        completed = run_lcfs_p(data_full, T_max);
    end
    
    % completed columns: [FinishTime, GenTime, UserID, MeasError]
    num_processed_total = size(completed, 1);
    
    % --- Inicializar saídas por usuário ---
    actuation_time_per_user = zeros(num_users, 1);
    processed_counts_per_user = zeros(num_users, 1);

    % --- Cálculo de Métricas se a fila não estiver vazia ---
    if isempty(completed)
        % Se nada foi processado, erro cresce linearmente
        v_mean = mean(velocities);
        avg_aeol_system = 0.5 * v_mean * T_max; 
        
        % Erro individual: Maximo
        avg_aeol_user = zeros(num_users, 1);
        for u = 1:num_users
            avg_aeol_user(u) = 0.5 * velocities(u) * T_max;
        end
        return; % Retorna com zeros nas contagens e tempos
    end

    % 1. Contagem de Pacotes por Usuário
    user_ids = completed(:, 3);
    for u = 1:num_users
        processed_counts_per_user(u) = sum(user_ids == u);
    end

    % 2. Tempo de Atuação por Usuário (Baseado na "Posse" do Monitor)
    % Ordenar por tempo de chegada no monitor (Finish Time)
    data_sys = sortrows(completed, 1); 
    finish_times = data_sys(:, 1);
    sorted_users = data_sys(:, 3);
    
    current_t = finish_times(1); % O sistema começa a ter info atualizada aqui
    
    % O tempo de 0 até o primeiro pacote não conta como "atuação" de nenhum usuário (é tempo "morto" ou info velha)
    
    for k = 1:length(finish_times)
        u_idx = sorted_users(k);
        
        % Determinar fim do intervalo deste pacote
        if k < length(finish_times)
            next_t = finish_times(k+1);
        else
            next_t = T_max; % O ultimo pacote reina até o fim da simulação
        end
        
        % Se o próximo evento passar de T_max (segurança), trava
        if next_t > T_max, next_t = T_max; end
        
        dt = next_t - current_t;
        
        if dt > 0
            actuation_time_per_user(u_idx) = actuation_time_per_user(u_idx) + dt;
            current_t = next_t;
        end
    end

    % 3. Cálculo AEoL do Sistema
    v_mean = mean(velocities);
    data_calc = completed(:, [1, 2, 4]); % [Finish, Gen, Error]
    avg_aeol_system = calc_position_error_area(data_calc, v_mean, T_max);

    % 4. Cálculo AEoL por Usuário
    avg_aeol_user = zeros(num_users, 1);
    for u = 1:num_users
        v_u = velocities(u);
        mask = (completed(:, 3) == u);
        data_u = completed(mask, [1, 2, 4]); 
        
        if isempty(data_u)
            avg_aeol_user(u) = 0.5 * v_u * T_max;
        else
            avg_aeol_user(u) = calc_position_error_area(data_u, v_u, T_max);
        end
    end
end

% =========================================================================
%  CORE FUNCTIONS (ÁREA, INTEGRAL, LCFS) - MANTIDAS IGUAIS
% =========================================================================
function avg_error = calc_position_error_area(data, v, T_max)
    data = sortrows(data, 1);
    T_recep = data(:, 1); T_gen = data(:, 2); Errors = data(:, 3);
    
    total_area = 0;
    t_current = 0; 
    last_gen_time = 0; last_error = 0; 
    
    for k = 1:length(T_recep)
        t_next = min(T_recep(k), T_max);
        
        if t_next > t_current
            a = v; b = v * last_gen_time + last_error;
            total_area = total_area + integrate_abs_linear(a, b, t_current, t_next);
            t_current = t_next;
        end
        if t_current >= T_max, break; end
        last_gen_time = T_gen(k); last_error = Errors(k);
    end
    
    if t_current < T_max
        a = v; b = v * last_gen_time + last_error;
        total_area = total_area + integrate_abs_linear(a, b, t_current, T_max);
    end
    avg_error = total_area / T_max;
end

function area = integrate_abs_linear(a, b, t1, t2)
    if t1 >= t2, area = 0; return; end
    t_z = b / a;
    if (t_z > t1) && (t_z < t2)
        area = integral_simple(a, b, t1, t_z) + integral_simple(a, b, t_z, t2);
    else
        area = integral_simple(a, b, t1, t2);
    end
end

function val = integral_simple(a, b, t_s, t_e)
    t_m = (t_s + t_e)/2; 
    if (a*t_m - b) < 0, s = -1; else, s = 1; end
    val = s * ((0.5*a*t_e^2 - b*t_e) - (0.5*a*t_s^2 - b*t_s));
end

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
        if fin_t <= T_max, completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)]; time_now = fin_t;
        else, time_now = fin_t; end
        served_mask(target) = true;
    end
end

function completed = run_lcfs_p(data, T_max)
    N = size(data, 1); completed = []; if N==0, return; end
    for k = 1:N-1
        ft = data(k,1) + data(k,3);
        if ft <= data(k+1,1) && ft <= T_max, completed = [completed; ft, data(k,1), data(k,2), data(k,4)]; end
    end
    ft = data(N,1) + data(N,3);
    if ft <= T_max, completed = [completed; ft, data(N,1), data(N,2), data(N,4)]; end
end