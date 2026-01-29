function [avg_aeol_user, avg_aeol_system, num_processed_total, processed_counts_per_user] = queue_logic_movement(arrivals, services, meas_errors, path_segments, type, T_max)
    
    % arrivals: [Time, UserID]
    % services: [ServiceTime]
    % meas_errors: [ErrorX, ErrorY]
    % path_segments: Trajetória real RWP
    
    num_users = max(arrivals(:,2));
    if isempty(num_users)
        num_users = 1;
    end
    
    data_full = [arrivals, services, meas_errors]; 

    % --- 1. Rodar Fila (Define quem é servido) ---
    if strcmp(type, 'NP')
        completed = run_lcfs_np(data_full, T_max);
    else
        completed = run_lcfs_p(data_full, T_max);
    end
    % completed: [FinishTime, GenTime, UserID, ErrX, ErrY]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filtra atualizações do Usuário 1 para visualizar
    user_id = 1;
    user_updates = completed(completed(:, 3) == user_id, :);

    % Chama a função de visualização
    debug_visualize_mse(user_updates, path_segments, T_max);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    num_processed_total = size(completed, 1);
    processed_counts_per_user = zeros(num_users, 1);
    
    if num_processed_total > 0
        for u = 1:num_users
            processed_counts_per_user(u) = sum(completed(:,3) == u);
        end
    end

    % --- PRE-CALCULO DO "ERRO CEGO" (BLIND ERROR) ---
    % Se não houver atualizações, o erro não é zero! 
    % O erro é a integral da divergência da trajetória em relação à posição inicial (p0).
    % Calculamos isso uma vez passando uma lista de updates vazia [].
    blind_mse = calculate_aeol([], path_segments, T_max);

    % --- TRATAMENTO DE CASO VAZIO TOTAL ---
    if isempty(completed)
        % Se NINGUÉM entregou pacote, todo mundo assume o erro cego
        avg_aeol_system = blind_mse;
        avg_aeol_user = ones(num_users, 1) * blind_mse;
        return; % Sai da função economizando tempo
    end
    
    % --- 2. Calcular AEoL do SISTEMA (União das Informações) ---
    % Usa TODOS os pacotes completados para estimar a posição
    avg_aeol_system = calculate_aeol(completed, path_segments, T_max);
    
    % --- 3. Calcular AEoL por USUÁRIO (Cálculo Separado) ---
    % Calcula o erro como se apenas as atualizações do Usuário U existissem
    avg_aeol_user = zeros(num_users, 1);
    
    for u = 1:num_users
        % Filtra apenas atualizações deste usuário
        user_updates = completed(completed(:, 3) == u, :);
        avg_aeol_user(u) = calculate_aeol(user_updates, path_segments, T_max);
    end

end

%% --- SUB-FUNÇÃO: CÁLCULO DO MSE (ÁREA) ---
function mse = calculate_aeol(updates, segments, T_max)
    % updates: matriz [FinishTime, GenTime, UserID, ErrX, ErrY]
    % segments: trajetória real [tstart, tend, xstart,ystart, vx, vy]
    
    total_error_area = 0;
    current_time = 0;
    
    % Posição Inicial: 
    p0 = get_position_at_t(0, segments);
    estimated_pos = p0;
    
    % Ordenar eventos de atualização
    if ~isempty(updates)
        events = sortrows(updates, 1); % Ordena por FinishTime
    else
        events = [];
    end
    
    event_idx = 1;
    num_events = size(events, 1);
    
    while current_time < T_max
        % Próximo evento ou fim da simulação
        if event_idx <= num_events
            next_time = min(events(event_idx, 1), T_max);
            is_update = true;
        else
            next_time = T_max;
            is_update = false;
        end
        
        % Integrar erro no intervalo [current, next]
        if next_time > current_time
            %segment_error = integrate_sq_error(current_time, next_time, estimated_pos, segments);
            segment_error = integrate_distance_error(current_time, next_time, estimated_pos, segments);
            total_error_area = total_error_area + segment_error;
        end
        
        current_time = next_time;
        
        % Atualizar Estimativa
        if is_update && current_time < T_max
            % Pega dados do pacote
            gen_time = events(event_idx, 2);
            err_x = events(event_idx, 4);
            err_y = events(event_idx, 5);
            
            % Nova estimativa = Posição Real na Geração + Erro do Sensor
            pos_at_gen = get_position_at_t(gen_time, segments);
            estimated_pos = pos_at_gen + [err_x, err_y];
            
            event_idx = event_idx + 1;
        end
    end
    
    mse = total_error_area / T_max;
end

%% --- Funções Auxiliares Geométricas ---

function pos = get_position_at_t(t, segments)
    idx = find(segments(:,1) <= t & segments(:,2) >= t, 1, 'last');
    if isempty(idx)
        idx = size(segments, 1); 
        t = min(t, segments(idx, 2));
    end
    row = segments(idx, :);
    dt = t - row(1);
    pos = [row(3) + row(5)*dt, row(4) + row(6)*dt];
end


function area = integrate_distance_error(t1, t2, P_est, segments)
    % Integra || P(t) - P_est || (Distância Euclidiana Linear)
    % Resultado em [Metros * Segundos]
    
    area = 0;
    t_curr = t1;
    
    while t_curr < t2
        % 1. Achar segmento atual da trajetória
        idx = find(segments(:,1) <= t_curr & segments(:,2) > t_curr, 1, 'last');
        if isempty(idx), idx = find(segments(:,2) >= t_curr, 1, 'first'); end
        if isempty(idx), break; end
        
        row = segments(idx, :);
        t_limit = min(t2, row(2));
        
        % Limites de integração relativos ao início do segmento
        tau_a = t_curr - row(1);
        tau_b = t_limit - row(1);
        
        % Dados do vetor movimento
        A = [row(3), row(4)]; % Posição inicial (x,y)
        V = [row(5), row(6)]; % Velocidade (vx,vy)
        C = A - P_est;        % Diferença inicial (offset)
        
        % 2. Definir a função da distância no tempo (handle function)
        % D(tau) = sqrt( (Cx + Vx*tau)^2 + (Cy + Vy*tau)^2 )
        dist_func = @(tau) sqrt( (C(1) + V(1).*tau).^2 + (C(2) + V(2).*tau).^2 );
        
        % 3. Integração Numérica Eficiente
        % O MATLAB calcula a área sob a curva da distância
        segment_area = integral(dist_func, tau_a, tau_b);
        
        area = area + segment_area;
        t_curr = t_limit;
    end
end

function area = integrate_sq_error(t1, t2, P_est, segments)
    % Integra || P(t) - P_est ||^2 analiticamente
    area = 0;
    t_curr = t1;
    
    while t_curr < t2
        % Achar segmento atual da trajetória
        idx = find(segments(:,1) <= t_curr & segments(:,2) > t_curr, 1, 'last');
        if isempty(idx), idx = find(segments(:,2) >= t_curr, 1, 'first'); end
        if isempty(idx), break; end
        
        row = segments(idx, :);
        t_limit = min(t2, row(2));
        
        % Integração Analítica
        tau_a = t_curr - row(1);
        tau_b = t_limit - row(1);
        
        A = [row(3), row(4)]; % Posição inicial do segmento
        V = [row(5), row(6)]; % Velocidade
        C = A - P_est;        % Diferença base
        
        % Coeficientes do polinômio quadrático de erro
        K0 = C(1)^2 + C(2)^2;
        K1 = 2 * (C(1)*V(1) + C(2)*V(2));
        K2 = V(1)^2 + V(2)^2;
        
        val_b = K0*tau_b + (K1/2)*tau_b^2 + (K2/3)*tau_b^3;
        val_a = K0*tau_a + (K1/2)*tau_a^2 + (K2/3)*tau_a^3;
        
        area = area + (val_b - val_a);
        t_curr = t_limit;
    end
end

%% --- Lógica de Fila LCFS ---

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
            time_now = data(upcoming, 1); continue;
        end
        target = queue(end);
        if length(queue) > 1, served_mask(queue(1:end-1)) = true; end
        
        fin_t = time_now + data(target, 3);
        if fin_t <= T_max
            completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4), data(target, 5)];
            time_now = fin_t;
        else
            time_now = T_max;
        end
        served_mask(target) = true;
    end
end

function completed = run_lcfs_p(data, T_max)
    N = size(data, 1);
    completed = [];
    if N == 0, return; end
    
    for k = 1:N-1
        ft = data(k,1) + data(k,3);
        % Preempção se o próximo chega antes deste terminar
        if ft <= data(k+1,1) && ft <= T_max 
            completed = [completed; ft, data(k,1), data(k,2), data(k,4), data(k,5)];
        end
    end
    k = N; 
    ft = data(k,1) + data(k,3);
    if ft <= T_max
        completed = [completed; ft, data(k,1), data(k,2), data(k,4), data(k,5)];
    end
end
