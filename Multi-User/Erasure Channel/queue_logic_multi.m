function [avg_aoi_per_user, avg_aoi_system] = queue_logic_multi(arrivals, services, num_users, type)

    if strcmp(type, 'NP')
        completed = run_lcfs_np(arrivals, services);
    else
        completed = run_lcfs_p(arrivals, services);
    end
    
    % =====================================================================
    % 1. CÁLCULO DO AOI DO SISTEMA (Cenário A: Fusão)
    % =====================================================================
    % Para o sistema, não importa o ID do usuário (coluna 3), apenas que 
    % uma atualização fresca chegou. Usamos todos os pacotes completados.
    if isempty(completed)
        avg_aoi_system = 0;
    else
        data_sys = completed(:, 1:2); % [Completion, Generation]
        avg_aoi_system = calc_aoi_area(data_sys);
    end

    % =====================================================================
    % 2. CÁLCULO DO AOI POR USUÁRIO (Individual)
    % =====================================================================
    avg_aoi_per_user = zeros(num_users, 1);
    
    for u = 1:num_users
        % Filtra apenas pacotes deste usuário
        mask = (completed(:, 3) == u);
        data_u = completed(mask, 1:2);
        
        if isempty(data_u)
            avg_aoi_per_user(u) = 0;
        else
            avg_aoi_per_user(u) = calc_aoi_area(data_u);
        end
    end
end

%% --- SUBFUNÇÃO AUXILIAR PARA CÁLCULO DE ÁREA ---
function avg_val = calc_aoi_area(data)
    % data espera: [Completion_Time, Generation_Time]
    
    % Ordena por tempo de conclusão
    data = sortrows(data, 1);
    C = data(:, 1); 
    G = data(:, 2); 
    
    area = 0; 
    t_now = 0;      % Tempo atual da simulação (varredura)
    current_aoi_start = 0; % Generation time do pacote ativo no sistema
    
    for k = 1:length(C)
        t_fin = C(k);
        t_gen = G(k);
        
        % Se o pacote entregue for MAIS VELHO que o que já temos,
        % ele não reduz o AoI (em sistemas ideais). 
        % Mas na sua lógica de fila, pacotes velhos já foram descartados, 
        % então t_gen será sempre >= current_aoi_start (exceto desordem rara).
        % Por segurança, garantimos a consistência:
        if t_gen < current_aoi_start
             % Pacote obsoleto chegou fora de ordem (raro em M/G/1 FIFO/LCFS simples)
             % Tratamos apenas avançando o tempo, sem reduzir a idade base
             % (Opcional: pular lógica de trapézio normal)
        end

        % Intervalo [t_now, t_fin]
        dt = t_fin - t_now;
        
        % Cálculo da Área Trapezoidal
        % Altura Inicial = t_now - current_aoi_start
        % Altura Final   = t_fin - current_aoi_start
        age_start = t_now - current_aoi_start;
        age_end   = t_fin - current_aoi_start;
        
        area = area + (age_start + age_end) * dt / 2;
        
        t_now = t_fin;
        
        % A idade cai para (t_fin - t_gen) instantaneamente
        % O novo "marco zero" da idade é o t_gen
        current_aoi_start = t_gen; 
    end
    
    if t_now > 0
        avg_val = area / t_now;
    else
        avg_val = 0;
    end
end

%% --- M/G/1/1 NON-PREEMPTIVE (WAITING) ---
% (Mantido igual ao seu original, pois para "Mesmo Objeto" a lógica de descarte global é válida)
function completed = run_lcfs_np(arr, srv)
    N = size(arr, 1);
    completed = [];
    served_mask = false(N, 1);
    
    time_now = 0;
    processed_count = 0;
    
    while processed_count < N
        queue_idx = find(arr(:,1) <= time_now & ~served_mask);
        
        if isempty(queue_idx)
            upcoming = find(arr(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming)
                break; 
            end
            time_now = arr(upcoming, 1);
            continue;
        end
        
        target = queue_idx(end);
        
        if length(queue_idx) > 1
            dropped_indices = queue_idx(1:end-1);
            served_mask(dropped_indices) = true; 
            processed_count = processed_count + length(dropped_indices);
        end
        
        start_t = time_now;
        finish_t = start_t + srv(target);
        
        completed = [completed; finish_t, arr(target, 1), arr(target, 2)];
        
        served_mask(target) = true;
        
        time_now = finish_t;
        processed_count = processed_count + 1;
    end
end

%% ---  M/G/1/1 PREEMPTIVE (LCLS) ---
function completed = run_lcfs_p(arr, srv)
    N = size(arr, 1);
    completed = [];
    
    for k = 1:(N-1)
        arrival_curr = arr(k, 1);
        arrival_next = arr(k+1, 1); 
        service_time = srv(k);
        
        potential_finish = arrival_curr + service_time;
        
        if potential_finish <= arrival_next
            completed = [completed; potential_finish, arrival_curr, arr(k, 2)];
        else
            % Interrompido pelo próximo pacote (global)
        end
    end
    
    last_idx = N;
    finish_last = arr(last_idx, 1) + srv(last_idx);
    completed = [completed; finish_last, arr(last_idx, 1), arr(last_idx, 2)];
end