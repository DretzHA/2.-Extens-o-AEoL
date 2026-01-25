function [avg_aeol_user, avg_aeol_system, num_processed_total, processed_counts_per_user] = queue_logic_trajectory(arrivals, services, errors, velocities, type, T_max)
    
    num_users = length(velocities);
    data_full = [arrivals, services, errors];

    % Rodar Fila
    if strcmp(type, 'NP')
        completed = run_lcfs_np(data_full, T_max);
    else
        completed = run_lcfs_p(data_full, T_max);
    end
    
    % completed: [FinishTime, GenTime, UserID, MeasError]
    num_processed_total = size(completed, 1);
    
    % --- Inicializar resultados por usuário ---
    processed_counts_per_user = zeros(num_users, 1);

    if isempty(completed)
        % Se nada foi processado, erro cresce linearmente
        v_mean = mean(velocities);
        avg_aeol_system = 0.5 * v_mean * T_max; 
        
        % Erro individual
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

    % 2. Cálculo AEoL do Sistema
    v_mean = mean(velocities);
    data_calc = completed(:, [1, 2, 4]); % [Finish, Gen, Error]
    avg_aeol_system = calc_position_error_area(data_calc, v_mean, T_max);

    % 3. Cálculo AEoL por Usuário
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
%  FUNCIO
% =========================================================================
function avg_error = calc_position_error_area(data, v, T_max)
    data = sortrows(data, 1);
    T_recep = data(:, 1); 
    T_gen = data(:, 2);
    Errors = data(:, 3);
    
    total_area = 0;
    t_current = 0; 
    last_gen_time = 0; last_error = 0; 
    
    for k = 1:length(T_recep)
        t_next = min(T_recep(k), T_max); %evita pegar tempo maior que o de simulação
        
        if t_next > t_current
            %a - inclinação da reta
            %b - estimativa atual
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

function area = integrate_abs_linear(a, b, t_current, t_next)
    if t_current >= t_next
        area = 0; 
        return; 
    end

    t_z = b / a; %ponto de corte (para caso com RMSE)

    if (t_z > t_current) && (t_z < t_next) % se há corte, faz calculo de duas areas
        area = integral_simple(a, b, t_current, t_z) + integral_simple(a, b, t_z, t_next);
    else
        area = integral_simple(a, b, t_current, t_next);
    end
end

function val = integral_simple(a, b, t_current, t_next)
    %calculo da area por integral
    t_m = (t_current + t_next)/2; 

    if (a*t_m - b) < 0
        s = -1; 
    else
        s = 1;
    end

    val = s * ((0.5*a*t_next^2 - b*t_next) - (0.5*a*t_current^2 - b*t_current));
end


%FUNÇÃO FILA NÃO PREEMPTIVA
function completed = run_lcfs_np(data, T_max)
    N = size(data, 1); 
    completed = []; %armazena dados completos
    served_mask = false(N, 1); %flag para verificar se já foi processado
    time_now = 0; 
    while true
        if time_now >= T_max %se ultrapassa tempo total, então para simulação
            break;
        end

        queue = find(data(:,1) <= time_now & ~served_mask); %procura primeiro valor da fila que não foi servido

        if isempty(queue)
            upcoming = find(data(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming) %se já foi todos
                break; 
            end
            time_now = data(upcoming, 1); 
            continue;
        end

        target = queue(end);

        if length(queue) > 1 %se existe mais uma informaçao na fila, somente deixa a mais atual (buffer 1)
            served_mask(queue(1:end-1)) = true;
        end

        fin_t = time_now + data(target, 3); %tempo de termino

        if fin_t <= T_max
            completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)]; 
            time_now = fin_t;
        else
            time_now = fin_t; 
        end
        served_mask(target) = true;
    end
end

function completed = run_lcfs_p(data, T_max)
    N = size(data, 1); 
    completed = []; 

    if N==0
        return; 
    end

    for k = 1:N-1
        ft = data(k,1) + data(k,3);

        if ft <= data(k+1,1) && ft <= T_max %verifica se atualização termina antes do proximo arrival
            completed = [completed; ft, data(k,1), data(k,2), data(k,4)];
        end

    end

    ft = data(N,1) + data(N,3); %ultima informação verifica pelo tempo máximo da simulação

    if ft <= T_max
        completed = [completed; ft, data(N,1), data(N,2), data(N,4)];
    end
end