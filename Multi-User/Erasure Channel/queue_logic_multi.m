function avg_aoi_per_user = queue_logic_multi(arrivals, services, num_users, type)

    if strcmp(type, 'NP')
        completed = run_lcfs_np(arrivals, services);
    else
        completed = run_lcfs_p(arrivals, services);
    end
    
    % AOI GRAPHICAL AREA BY USER
    avg_aoi_per_user = zeros(num_users, 1);
    
    for u = 1:num_users
        % Obtained the completed updates by each user
        mask = (completed(:, 3) == u);
        data = completed(mask, 1:2);
        
        if isempty(data)
            avg_aoi_per_user(u) = 0;
            continue;
        end
        
        % Sort by time
        data = sortrows(data, 1);
        C = data(:, 1); % Completion times
        G = data(:, 2); % Generation times
        
        area = 0; 
        t_prev = 0; % Initial Aoi equal to zero
        last_gen = 0;
           
        
        % Calculo por Area Trapezio
        current_aoi_start = 0;
        t_now = 0;
        
        for k = 1:length(C)
            t_fin = C(k);
            t_gen = G(k);
            
            % Interval [t_now, t_fin]
            dt = t_fin - t_now;
            
            % Area do trapézio: Base_menor + Base_maior * Altura / 2
            % Base menor = t_now - current_aoi_start
            % Base maior = t_fin - current_aoi_start
            age_start = t_now - current_aoi_start;
            age_end   = t_fin - current_aoi_start;
            
            area = area + (age_start + age_end) * dt / 2;
            
            t_now = t_fin;
            current_aoi_start = t_gen; 
        end
        
        avg_aoi_per_user(u) = area / t_now;
    end
end

%% --- M/G/1/1 NON-PREEMPTIVE (WAITING) ---
function completed = run_lcfs_np(arr, srv)
    N = size(arr, 1);
    completed = [];
    served_mask = false(N, 1);
    
    time_now = 0;
    processed_count = 0;
    
    while processed_count < N
        % Encontra todos os pacotes que chegaram até agora e não foram servidos
        % nem descartados
        queue_idx = find(arr(:,1) <= time_now & ~served_mask);
        
        if isempty(queue_idx)
            % Avança o tempo para a próxima chegada se a fila estiver vazia
            upcoming = find(arr(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming)
                break; % Acabaram os pacotes
            end
            time_now = arr(upcoming, 1);
            continue;
        end
        
        % LÓGICA ALTERADA: M/G/1/1 com Descarte de Obsoletos
        % Seleciona o pacote MAIS RECENTE (último índice)
        target = queue_idx(end);
        
        % DESCARTA todos os pacotes anteriores que estavam na fila
        % (Ignora informações anteriores à última atualização)
        if length(queue_idx) > 1
            dropped_indices = queue_idx(1:end-1);
            served_mask(dropped_indices) = true; % Marca como "processado" (ignorado)
            processed_count = processed_count + length(dropped_indices);
        end
        
        % Processa o alvo
        start_t = time_now;
        finish_t = start_t + srv(target);
        
        % Registra conclusão: [Finish, Generation, UserID]
        completed = [completed; finish_t, arr(target, 1), arr(target, 2)];
        
        served_mask(target) = true;
        
        % O servidor fica ocupado até finish_t
        time_now = finish_t;
        processed_count = processed_count + 1;
    end
end

%% ---  M/G/1/1 PREEMPTIVE (LCLS) ---
function completed = run_lcfs_p(arr, srv)
    N = size(arr, 1);
    completed = [];
    
    % Na lógica preemptiva pura (M/G/1/1 Preemptive), qualquer chegada nova
    % interrompe (mata) o serviço atual imediatamente.
    % O código abaixo verifica se o pacote 'k' consegue terminar antes de 'k+1' chegar.
    
    for k = 1:(N-1)
        arrival_curr = arr(k, 1);
        arrival_next = arr(k+1, 1); % Próxima chegada (interrupção potencial)
        service_time = srv(k);
        
        potential_finish = arrival_curr + service_time;
        
        if potential_finish <= arrival_next
            % Sucesso: Terminou antes de ser interrompido
            completed = [completed; potential_finish, arrival_curr, arr(k, 2)];
        else
            % Falha: Foi interrompido. O pacote 'k' é descartado.
            % O sistema passa a servir o 'k+1' imediatamente no tempo arr(k+1,1).
            % (Não fazemos nada aqui, apenas não adicionamos ao 'completed')
        end
    end
    
    % O último pacote sempre completa (ninguém para interromper)
    last_idx = N;
    finish_last = arr(last_idx, 1) + srv(last_idx);
    completed = [completed; finish_last, arr(last_idx, 1), arr(last_idx, 2)];
end