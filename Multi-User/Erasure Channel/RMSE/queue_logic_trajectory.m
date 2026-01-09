function [avg_aeol_user, avg_aeol_system, num_processed] = queue_logic_trajectory(arrivals, services, errors, velocities, type, T_max)
    % arrivals: [Time, UserID]
    % services: [Duration]
    % errors:   [MeasurementError]
    % T_max:    Tempo total da simulação
    
    data_full = [arrivals, services, errors];

    % Rodar Fila
    if strcmp(type, 'NP')
        completed = run_lcfs_np(data_full, T_max);
    else
        completed = run_lcfs_p(data_full, T_max);
    end
    
    % completed: [FinishTime, GenTime, UserID, MeasError]
    num_processed = size(completed, 1);

    % --- Cálculo AEoL do Sistema ---
    if isempty(completed)
        % Se nada foi processado, o erro cresce linearmente de 0 a T_max
        % Integral de |v*t| dt = 0.5 * v * T^2 -> Média = 0.5 * v * T
        v_mean = mean(velocities);
        avg_aeol_system = 0.5 * v_mean * T_max; 
    else
        v_mean = mean(velocities);
        data_sys = completed(:, [1, 2, 4]); % [Finish, Gen, Error]
        avg_aeol_system = calc_position_error_area(data_sys, v_mean, T_max);
    end

    % --- Cálculo AEoL por Usuário (CORREÇÃO AQUI) ---
    avg_aeol_user = zeros(length(velocities), 1);
    
    if isempty(completed)
        % CASO 1: Matriz vazia. Nenhum pacote processado para NINGUÉM.
        % Atribui o erro máximo (crescimento linear) para todos.
        for u = 1:length(velocities)
            avg_aeol_user(u) = 0.5 * velocities(u) * T_max;
        end
    else
        % CASO 2: Existem pacotes. Filtramos por usuário.
        for u = 1:length(velocities)
            v_u = velocities(u);
            
            % Agora é seguro acessar a coluna 3
            mask = (completed(:, 3) == u);
            data_u = completed(mask, [1, 2, 4]); % [Finish, Gen, Error]
            
            if isempty(data_u)
                % Usuário específico não teve pacotes, mas outros tiveram
                avg_aeol_user(u) = 0.5 * v_u * T_max;
            else
                avg_aeol_user(u) = calc_position_error_area(data_u, v_u, T_max);
            end
        end
    end
end

% =========================================================================
%  CORE DO CÁLCULO DE ÁREA (INTEGRAL DO ERRO DE TRAJETÓRIA)
% =========================================================================
function avg_error = calc_position_error_area(data, v, T_max)
    % data: [ReceptionTime, GenerationTime, MeasError]
    % v: Velocidade do alvo
    % T_max: Limite superior da integral
    
    data = sortrows(data, 1);
    
    T_recep = data(:, 1);
    T_gen   = data(:, 2);
    Errors  = data(:, 3);
    
    total_area = 0;
    t_current = 0; 
    
    last_gen_time = 0;
    last_error    = 0; 
    
    % Integra intervalor entre pacotes
    for k = 1:length(T_recep)
        t_next_event = T_recep(k);
        
        % Se o próximo evento ocorre APÓS o fim da simulação (segurança), trava.
        if t_next_event > T_max
            t_next_event = T_max;
        end
        
        if t_next_event > t_current
            a = v;
            b = v * last_gen_time + last_error;
            seg_area = integrate_abs_linear(a, b, t_current, t_next_event);
            total_area = total_area + seg_area;
            t_current = t_next_event;
        end
        
        % Se já alcançamos o fim da simulação, pare.
        if t_current >= T_max
            break;
        end

        % Atualiza estado
        last_gen_time = T_gen(k);
        last_error    = Errors(k);
    end
    
    % Integra o "Resto" (do último pacote até T_max)
    if t_current < T_max
        a = v;
        b = v * last_gen_time + last_error;
        seg_area = integrate_abs_linear(a, b, t_current, T_max);
        total_area = total_area + seg_area;
    end
    
    avg_error = total_area / T_max;
end

function area = integrate_abs_linear(a, b, t1, t2)
    if t1 >= t2, area = 0; return; end
    t_zero = b / a;
    if (t_zero > t1) && (t_zero < t2)
        area = integral_simple(a, b, t1, t_zero) + integral_simple(a, b, t_zero, t2);
    else
        area = integral_simple(a, b, t1, t2);
    end
end

function val = integral_simple(a, b, t_start, t_end)
    t_mid = (t_start + t_end) / 2;
    f_mid = a * t_mid - b;
    F_end   = 0.5 * a * t_end^2   - b * t_end;
    F_start = 0.5 * a * t_start^2 - b * t_start;
    raw_integral = F_end - F_start;
    if f_mid < 0, val = -raw_integral; else, val = raw_integral; end
end

% --- LCFS-NP (Não-Preemptivo) ---
function completed = run_lcfs_np(data, T_max)
    N = size(data, 1);
    completed = [];
    served_mask = false(N, 1);
    time_now = 0;
    
    while true
        if time_now >= T_max, break; end
        
        queue = find(data(:,1) <= time_now & ~served_mask);
        
        if isempty(queue)
            upcoming = find(data(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming), break; end
            time_now = data(upcoming, 1);
            continue;
        end
        
        target = queue(end);
        
        if length(queue) > 1
            drp = queue(1:end-1);
            served_mask(drp) = true;
        end
        
        fin_t = time_now + data(target, 3);
        
        if fin_t <= T_max
            completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)];
            time_now = fin_t;
        else
            time_now = fin_t; 
        end
        served_mask(target) = true;
    end
end

% --- LCFS-P (Preemptivo) ---
function completed = run_lcfs_p(data, T_max)
    N = size(data, 1);
    completed = [];
    
    if N == 0, return; end
    
    for k = 1:N-1
        arr_c = data(k,1);
        arr_n = data(k+1,1);
        serv  = data(k,3);
        
        if arr_c + serv <= arr_n
            fin_t = arr_c + serv;
            if fin_t <= T_max
                completed = [completed; fin_t, arr_c, data(k,2), data(k,4)];
            end
        end
    end
    
    last = N;
    arr_last = data(last, 1);
    serv_last = data(last, 3);
    fin_last = arr_last + serv_last;
    
    if fin_last <= T_max
        completed = [completed; fin_last, arr_last, data(last, 2), data(last, 4)];
    end
end