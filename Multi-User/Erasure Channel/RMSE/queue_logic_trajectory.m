function [avg_aeol_user, avg_aeol_system] = queue_logic_trajectory(arrivals, services, errors, velocities, type)
    % arrivals: [Time, UserID]
    % services: [Duration]
    % errors:   [MeasurementError] (Mesma ordem de arrivals)
    % velocities: Vetor de velocidades por usuário
    
    % Montar matriz completa: [ArrTime, UserID, ServTime, Error]
    data_full = [arrivals, services, errors];

    % Rodar Fila (Lógica de aceitação/descarte não muda)
    if strcmp(type, 'NP')
        completed = run_lcfs_np(data_full);
    else
        completed = run_lcfs_p(data_full);
    end
    
    % completed: [FinishTime, GenTime, UserID, MeasError]

    % --- Cálculo AEoL do Sistema ---
    if isempty(completed)
        avg_aeol_system = 0;
    else
        % Para o sistema, consideramos a média dos erros dos usuários ativos
        % Ou, simplificando, tratamos como um fluxo único com velocidade média
        % (Depende da definição exata de "System AEoL". Aqui uso média ponderada).
        v_mean = mean(velocities);
        
        % Vamos calcular o erro acumulado total da trajetória "combinada"
        data_sys = completed(:, [1, 2, 4]); % [Finish, Gen, Error]
        avg_aeol_system = calc_position_error_area(data_sys, v_mean);
    end

    % --- Cálculo AEoL por Usuário ---
    avg_aeol_user = zeros(length(velocities), 1);
    for u = 1:length(velocities)
        v_u = velocities(u);
        mask = (completed(:, 3) == u);
        data_u = completed(mask, [1, 2, 4]); % [Finish, Gen, Error]
        
        if isempty(data_u)
            avg_aeol_user(u) = 0;
        else
            avg_aeol_user(u) = calc_position_error_area(data_u, v_u);
        end
    end
end

% =========================================================================
%  CORE DO CÁLCULO DE ÁREA (INTEGRAL DO ERRO DE TRAJETÓRIA)
% =========================================================================
function avg_error = calc_position_error_area(data, v)
    % data: [ReceptionTime, GenerationTime, MeasError]
    % v: Velocidade do alvo
    
    % Ordenar pelos tempos de recepção (mudança de estado no monitor)
    data = sortrows(data, 1);
    
    T_recep = data(:, 1);
    T_gen   = data(:, 2);
    Errors  = data(:, 3);
    
    total_area = 0;
    t_current = 0; % Tempo atual da simulação (início em 0)
    
    % Estado inicial do monitor (antes do primeiro pacote)
    % Assumimos erro 0 ou estado inicial conhecido? 
    % Normalmente em AoI, assume-se custo zero ou linear partindo de 0.
    % Vamos assumir que começamos em t=0 na posição 0 com estimativa 0.
    last_gen_time = 0;
    last_error    = 0; 
    
    for k = 1:length(T_recep)
        t_next_event = T_recep(k);
        
        % Intervalo [t_current, t_next_event]
        % Neste intervalo, o monitor mantém a informação do pacote ANTERIOR.
        % Posição Real: P(t) = v * t
        % Posição Est:  P_hat(t) = v * last_gen_time + last_error
        % Erro(t) = | P(t) - P_hat(t) | 
        %         = | v*t - (v*last_gen_time + last_error) |
        %         = | v*(t - last_gen_time) - last_error |
        
        % Definindo constantes para a integral |at - b|
        % a = v
        % b = v * last_gen_time + last_error
        % Queremos integral de |a*t - b| de t_current até t_next_event
        
        a = v;
        b = v * last_gen_time + last_error;
        
        seg_area = integrate_abs_linear(a, b, t_current, t_next_event);
        total_area = total_area + seg_area;
        
        % Atualiza estado para o próximo intervalo
        t_current = t_next_event;
        last_gen_time = T_gen(k);
        last_error    = Errors(k);
    end
    
    if t_current > 0
        avg_error = total_area / t_current;
    else
        avg_error = 0;
    end
end

function area = integrate_abs_linear(a, b, t1, t2)
    % Calcula integral de |a*t - b| dt de t1 a t2
    % A função f(t) = a*t - b é uma reta.
    % A raiz é t_zero = b/a.
    
    if t1 >= t2
        area = 0; return;
    end
    
    t_zero = b / a;
    
    % Verifica se a reta cruza o zero DENTRO do intervalo
    if (t_zero > t1) && (t_zero < t2)
        % Divide em duas partes: [t1, t_zero] e [t_zero, t2]
        area = integral_simple(a, b, t1, t_zero) + integral_simple(a, b, t_zero, t2);
    else
        % Intervalo direto
        area = integral_simple(a, b, t1, t2);
    end
end

function val = integral_simple(a, b, t_start, t_end)
    % Integral definida de |at - b| sabendo que não cruza zero (sinal constante)
    % Primitiva de (at - b) é 0.5*a*t^2 - b*t
    % Se a função for negativa no intervalo, invertemos o sinal do resultado.
    
    % Avaliar no ponto médio para saber o sinal
    t_mid = (t_start + t_end) / 2;
    f_mid = a * t_mid - b;
    
    F_end   = 0.5 * a * t_end^2   - b * t_end;
    F_start = 0.5 * a * t_start^2 - b * t_start;
    
    raw_integral = F_end - F_start;
    
    if f_mid < 0
        val = -raw_integral;
    else
        val = raw_integral;
    end
end

% --- Copiar as funções de fila LCFS (iguais ao anterior, só repassando erro) ---
function completed = run_lcfs_np(data)
    % data: [Arr, User, Serv, Err]
    N = size(data, 1);
    completed = [];
    served_mask = false(N, 1);
    time_now = 0;
    processed = 0;
    while processed < N
        queue = find(data(:,1) <= time_now & ~served_mask);
        if isempty(queue)
            upcoming = find(data(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming), break; end
            time_now = data(upcoming, 1);
            continue;
        end
        target = queue(end); % LCFS
        % Drop others
        if length(queue) > 1
            drp = queue(1:end-1);
            served_mask(drp) = true;
            processed = processed + length(drp);
        end
        
        fin_t = time_now + data(target, 3);
        completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)];
        served_mask(target) = true;
        time_now = fin_t;
        processed = processed + 1;
    end
end

function completed = run_lcfs_p(data)
    N = size(data, 1);
    completed = [];
    for k = 1:N-1
        arr_c = data(k,1); arr_n = data(k+1,1); serv = data(k,3);
        if arr_c + serv <= arr_n
            completed = [completed; arr_c+serv, arr_c, data(k,2), data(k,4)];
        end
    end
    last = N;
    completed = [completed; data(last,1)+data(last,3), data(last,1), data(last,2), data(last,4)];
end