function avg_aoi_per_user = queue_logic_multi(arrivals, services, num_users, type)
% QUEUE_LOGIC_MULTI Simula a disputa na fila compartilhada.
% Entradas:
%   arrivals: Matriz [Tempo, UserID] ordenada
%   services: Vetor de tempos de serviço
%   type: 'NP' (Sem Preempção) ou 'P' (Com Preempção)

    % Executa a simulação da fila dependendo do tipo
    if strcmp(type, 'NP')
        completed = run_lcfs_np(arrivals, services);
    else
        completed = run_lcfs_p(arrivals, services);
    end
    
    % CÁLCULO DA ÁREA SOB A CURVA (AoI) PARA CADA UTILIZADOR
    avg_aoi_per_user = zeros(num_users, 1);
    
    for u = 1:num_users
        % Extrai apenas os pacotes concluídos do utilizador 'u'
        % Colunas: [Tempo_Final, Tempo_Geração, UserID]
        mask = (completed(:, 3) == u);
        data = completed(mask, 1:2);
        
        if isempty(data)
            avg_aoi_per_user(u) = 0; 
            continue;
        end
        
        % Garante ordenação por tempo de chegada no destino
        data = sortrows(data, 1);
        C = data(:, 1); % Completion times
        G = data(:, 2); % Generation times
        
        area = 0; 
        t_prev = 0; % Tempo da última atualização no monitor
        g_prev = 0; % Geração da última atualização
        max_g = 0;  % Para filtrar pacotes obsoletos
        
        for i = 1:length(C)
            % Só processa se a informação for mais fresca que a atual
            if G(i) > max_g
                interval = C(i) - t_prev;
                
                % Altura inicial do trapézio
                h1 = t_prev - g_prev;
                % Altura final do trapézio
                h2 = C(i) - g_prev;
                
                area = area + 0.5 * (h1 + h2) * interval;
                
                t_prev = C(i);
                g_prev = G(i);
                max_g = G(i);
            end
        end
        
        % Média temporal
        if t_prev > 0
            avg_aoi_per_user(u) = area / t_prev;
        end
    end
end

%% --- LÓGICA INTERNA: LCFS SEM PREEMPÇÃO ---
function completed = run_lcfs_np(arr, srv)
    N = size(arr, 1);
    completed = []; 
    
    time_now = 0;
    server_free = 0;
    served_mask = false(N, 1);
    processed_count = 0;
    
    while processed_count < N
        % Avança o tempo se o servidor estiver livre e ninguém chegou ainda
        if time_now < server_free
            time_now = server_free;
        end
        
        % Quem está na fila neste momento?
        queue_idx = find(arr(:,1) <= time_now & ~served_mask);
        
        if isempty(queue_idx)
            % Ninguém na fila, salta para a próxima chegada
            next_idx = find(~served_mask, 1, 'first');
            if isempty(next_idx), break; end
            time_now = arr(next_idx, 1);
        else
            % LCFS: Escolhe o último que chegou (mais recente)
            % Aqui está a disputa: pode ser de QUALQUER utilizador
            target = queue_idx(end);
            
            start_t = time_now;
            finish_t = start_t + srv(target);
            
            % Regista: [Finish, Generation, UserID]
            completed = [completed; finish_t, arr(target, 1), arr(target, 2)];
            
            served_mask(target) = true;
            server_free = finish_t;
            time_now = finish_t;
            processed_count = processed_count + 1;
        end
    end
end

%% --- LÓGICA INTERNA: LCFS COM PREEMPÇÃO ---
function completed = run_lcfs_p(arr, srv)
    N = size(arr, 1);
    completed = [];
    
    % Em sistema partilhado com preempção global:
    % Um pacote é servido SE terminar antes de QUALQUER outro pacote chegar.
    
    for k = 1:(N-1)
        arrival_curr = arr(k, 1);
        arrival_next = arr(k+1, 1); % Próxima chegada (pode ser outro user!)
        service_time = srv(k);
        
        potential_finish = arrival_curr + service_time;
        
        if potential_finish <= arrival_next
            % Sucesso: Terminou antes de ser interrompido
            completed = [completed; potential_finish, arrival_curr, arr(k, 2)];
        else
            % Falha: Foi interrompido pela chegada k+1
        end
    end
    
    % O último pacote da simulação assume-se que termina
    completed = [completed; arr(N,1)+srv(N), arr(N,1), arr(N,2)];
end