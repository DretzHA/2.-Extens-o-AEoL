function avg_aoi = calculate_aoi_lcfs(inter_arrival_times, service_times, type)
% CALCULATE_AOI_LCFS Calcula o AoI médio para fila LCFS.
%   type: 'Preemptive' ou 'NonPreemptive'

    generation_times = cumsum(inter_arrival_times);
    
    % Seleciona a lógica baseada no tipo
    switch type
        case 'NonPreemptive'
            [C, G] = simulate_non_preemptive(generation_times, service_times);
        case 'Preemptive'
            [C, G] = simulate_preemptive(generation_times, service_times);
        otherwise
            error('Tipo inválido.');
    end
    
    % Cálculo Comum da Área (AoI Médio)
    if isempty(C)
        avg_aoi = 0;
        return;
    end
    
    % Ordena por tempo de chegada no destino para cálculo da área
    results = sortrows([C, G], 1);
    C_sorted = results(:, 1);
    G_sorted = results(:, 2);
    
    area_total = 0;
    completion_time_prev = 0;
    generation_time_prev = 0;
    max_gen_seen = 0;
    
    for k = 1:length(C_sorted)
        C_curr = C_sorted(k);
        G_curr = G_sorted(k);
        
        % Só conta se trouxer informação nova
        if G_curr > max_gen_seen
            interval = C_curr - completion_time_prev;
            H1 = completion_time_prev - generation_time_prev;
            H2 = C_curr - generation_time_prev;
            
            area_step = 0.5 * (H1 + H2) * interval;
            area_total = area_total + area_step;
            
            completion_time_prev = C_curr;
            generation_time_prev = G_curr;
            max_gen_seen = G_curr;
        end
    end
    
    if completion_time_prev > 0
        avg_aoi = area_total / completion_time_prev;
    else
        avg_aoi = 0;
    end
end

%% === FUNÇÕES AUXILIARES LOCAIS ===

function [completed_times, gen_times_out] = simulate_non_preemptive(G, S)
    N = length(G);
    completed_times = zeros(N, 1);
    gen_times_out = G;
    served_flag = false(N, 1);
    current_time = 0;
    count_served = 0;
    
    while count_served < N
        % Quem está na fila?
        queue_indices = find(G <= current_time & ~served_flag);
        
        if isempty(queue_indices)
            next_idx = find(~served_flag, 1, 'first');
            if isempty(next_idx), break; end
            current_time = G(next_idx);
        else
            % LCFS: Pega o último
            target_idx = queue_indices(end);
            finish_time = current_time + S(target_idx);
            completed_times(target_idx) = finish_time;
            served_flag(target_idx) = true;
            current_time = finish_time;
            count_served = count_served + 1;
        end
    end
    % Remove zeros (pacotes não processados por algum erro de borda)
    mask = completed_times > 0;
    completed_times = completed_times(mask);
    gen_times_out = gen_times_out(mask);
end

function [completed_times, gen_times_out] = simulate_preemptive(G, S)
    N = length(G);
    completed_times = [];
    gen_times_out = [];
    
    % Lógica LCFS Com Preempção (Descarte)
    for k = 1:(N-1)
        potential_finish = G(k) + S(k);
        next_arrival = G(k+1);
        
        if potential_finish <= next_arrival
            completed_times = [completed_times; potential_finish];
            gen_times_out = [gen_times_out; G(k)];
        end
    end
    % Último pacote
    completed_times = [completed_times; G(N) + S(N)];
    gen_times_out = [gen_times_out; G(N)];
end