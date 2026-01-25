clc;
clear;
close all;

%% ========================================================================

%  =======================================================================

% --- Parâmetros ---
mu = 1; %taxa serviço
velocities = [1]; %velocidades
rho = 0.05:0.05:1.2; %taxa de utilização
lambda_vec = rho * mu; %taxa de chegada
N_updates = 1e2; %numero de atualizações
N_MC = 500; %realizações para monte carlo

% --- Vetores para resultados ---
AEoL_M_woP_Analytic = zeros(length(velocities), length(lambda_vec));
AEoL_M_woP_Sim = zeros(length(velocities), length(lambda_vec));

AEoL_M_withP_Analytic = zeros(length(velocities), length(lambda_vec));
AEoL_M_withP_Sim = zeros(length(velocities), length(lambda_vec));

for j = 1:length(lambda_vec) %para cada valor da taxa lambda
    lambda = lambda_vec(j); %valor atual lambda
    cur_rho = rho(j); %valor de atualização
    
    %Valor  Analítico
    for i = 1:length(velocities)
        v = velocities(i); %velocidade atual        
        % M/M/1 without Preemption
        num = cur_rho^3*(1+cur_rho)^2 + cur_rho^4 + (1+cur_rho)^3;
        den = cur_rho^4*(1+cur_rho) + cur_rho*(1+cur_rho)^3;
        termo_1 = 1 + num/den;
        AEoL_M_woP_Analytic(i, j) = (v / mu) * termo_1;

        % M/M/1 with Preemption and one user
        AEoL_M_withP_Analytic(i, j) = v * (1/lambda + 1/mu);
    end
    

   %SIMULAÇÃO MONTE CARLO

   % Simulação de Monte Carlo
    avg_AoI_woP_M_acc = 0;
    avg_AoI_withP_M_acc = 0;
    
    for k = 1:N_MC

        % Simulação M/M/1 - leva em consideranção distribuição exponencial
        % para aleatoriedade das taxas das filas, executado em Monte Carlo
        inter_T_M = exprnd(1/lambda, N_updates, 1); %tempo entre chegada
        ServiceTime_M = exprnd(1/mu, N_updates, 1); %tempo de serviço

        avg_AoI_woP_M_acc = avg_AoI_woP_M_acc + calculate_avg_aoi(inter_T_M, ServiceTime_M); %calculo do AoI médio sem preempção
        avg_AoI_withP_M_acc = avg_AoI_withP_M_acc + calculate_avg_aoi_preemptive(inter_T_M, ServiceTime_M); %calculo do AoI médio com preempção
        
        
    end

    mean_AoI_woP_M = avg_AoI_woP_M_acc / N_MC; 
    mean_AoI_withP_M = avg_AoI_withP_M_acc / N_MC; 
   
    for i = 1:length(velocities)
        v = velocities(i);
        AEoL_M_woP_Sim(i, j) = v * mean_AoI_woP_M;
        AEoL_M_withP_Sim(i, j) = v * mean_AoI_withP_M;
    end
end


%% ========================================================================
%  PLOTAGEM DOS RESULTADOS
%  ========================================================================
custom_colors = {'#D95319', '#7E2F8E', '#000000'};

% --- FIGURA: M/M/1 ---
figure();
hold on;
set(gca, 'FontSize', 12);

for i = 1:length(velocities)
    plot(lambda_vec, AEoL_M_woP_Analytic(i, :), '-', ...
        'LineWidth', 1.5, 'Color', custom_colors{1}, ...
        'DisplayName', sprintf('Analyt. woP v=%.1f', velocities(i)));

    plot(lambda_vec, AEoL_M_woP_Sim(i, :), '*', ...
        'LineWidth', 1.5, 'Color', custom_colors{1}, ...
        'DisplayName', sprintf('Sim. woP v=%.1f', velocities(i)));

    plot(lambda_vec, AEoL_M_withP_Analytic(i, :), '--', ...
        'LineWidth', 1.5, 'Color', custom_colors{2}, ...
        'DisplayName', sprintf('Analyt. wP v=%.1f', velocities(i)));

    plot(lambda_vec, AEoL_M_withP_Sim(i, :), 'o', ...
        'LineWidth', 1.5, 'Color', custom_colors{2}, ...
        'DisplayName', sprintf('Sim. wP v=%.1f', velocities(i)));

end


xlabel('Taxa de Chegada \lambda', 'FontSize', 16);
ylabel('AEoL [m]', 'FontSize', 16);
legend('show', 'Location', 'north', 'FontSize', 16);

grid on; box on;
hold off;


%% ========================================================================
%  FUNÇÕES AUXILIARES
%  ========================================================================
% Cálculo de AoI Médio
function avg_aoi = calculate_avg_aoi(inter_arrival_times, service_times)
    generation_times = cumsum(inter_arrival_times); % Tempo de geração
    N = length(generation_times); 
    completion_times = zeros(N, 1); % Vetor para armazenar quando cada pacote termina
    
    served_flag = false(N, 1); % Controle de quais pacotes já foram atendidos
    current_time = 0;          % Relógio do servidor
    count_served = 0;

    % --- 1. SIMULAÇÃO DA FILA LCFS SEM PREEMPÇÃO ---
    while count_served < N
        % Encontra pacotes que já chegaram (G <= tempo atual) e ainda não foram atendidos
        queue_indices = find(generation_times <= current_time & ~served_flag);
        
        if isempty(queue_indices)
            % SE A FILA ESTÁ VAZIA:
            % O servidor avança no tempo até a chegada do próximo pacote não atendido
            next_arrival_idx = find(~served_flag, 1, 'first');
            if isempty(next_arrival_idx)
                break; % Todos atendidos (segurança)
            end
            current_time = generation_times(next_arrival_idx);
            % Volta para o início do loop para processar este pacote
            continue; 
        else
            % SE TEM GENTE NA FILA:
            % LCFS: Escolhe o último que chegou (maior índice temporal)
            target_idx = queue_indices(end); 
            
            % Processa o pacote
            start_service = current_time;
            finish_service = start_service + service_times(target_idx);
            
            completion_times(target_idx) = finish_service;
            served_flag(target_idx) = true;
            
            % Servidor fica ocupado até terminar este serviço
            current_time = finish_service;
            count_served = count_served + 1;
        end
    end
    
    % --- 2. PREPARAÇÃO PARA CÁLCULO DO AOI (FILTRAGEM) ---
    % Como é LCFS, pacotes podem terminar fora de ordem.
    % O Monitor só atualiza se a informação for MAIS NOVA que a atual.
    
    % Tabela com [Tempo_Finalizacao, Tempo_Geracao]
    results = [completion_times, generation_times];
    
    % Ordena cronologicamente pelos tempos de finalização (chegada no destino)
    results = sortrows(results, 1);
    
    C_sorted = results(:, 1);
    G_sorted = results(:, 2);
    
    % --- 3. CÁLCULO DA ÁREA SOB A CURVA ---
    area_total = 0;
    completion_time_prev = 0;
    generation_time_prev = 0; % T0 assume-se geração 0 no tempo 0
    max_generation_seen = 0;  % Para garantir que não usaremos informação velha
    
    % O AoI começa crescendo linearmente desde t=0 até o primeiro pacote chegar
    % (Área inicial triangular ou trapezoidal dependendo da condição inicial)
    % Assumindo sistema vazio em t=0 com Age=0 ou Age inicial. 
    % O código original assume t_prev=0 e gen_prev=0.
    
    for k = 1:N
        C_curr = C_sorted(k);
        G_curr = G_sorted(k);
        
        % Só consideramos o pacote se ele trouxer informação mais nova (Fresher)
        if G_curr > max_generation_seen
            
            % Cálculo da área do trapézio entre a atualização anterior e a atual
            interval = C_curr - completion_time_prev;
            
            % Altura no início do intervalo (logo após a atualização anterior)
            H1 = completion_time_prev - generation_time_prev;
            
            % Altura no fim do intervalo (imediatamente antes da nova atualização)
            H2 = C_curr - generation_time_prev;
            
            area_step = 0.5 * (H1 + H2) * interval;
            area_total = area_total + area_step;
            
            % Atualiza variáveis para o próximo passo
            completion_time_prev = C_curr;
            generation_time_prev = G_curr; % O AoI cai para (C_curr - G_curr)
            max_generation_seen = G_curr;
        else
            % Se G_curr < max_generation_seen, o pacote é obsoleto (old news).
            % Ele chegou no destino, mas a informação já era velha comparada 
            % com o que o monitor já tinha. Não altera a curva do AoI (apenas 
            % continua subindo linearmente até o próximo pacote útil).
            % A área será computada no próximo pacote válido que englobará esse tempo.
        end
    end
    
    % Adiciona a área final residual (do último pacote até o fim da simulação se necessário)
    % O código original divide pelo completion_times(end). Vamos manter o padrão:
    % A área calculada até agora vai exatamente até o último completion_time útil.
    
    avg_aoi = area_total / completion_time_prev; 
end


%% ========================================================================
%  FUNÇÃO ALTERADA: LCFS COM PREEMPÇÃO
%  ========================================================================
function avg_aoi = calculate_avg_aoi_preemptive(inter_arrival_times, service_times)
    generation_times = cumsum(inter_arrival_times); % Tempo de geração (G)
    N = length(generation_times); 
    
    % Vetores para armazenar apenas os pacotes que tiveram sucesso
    successful_completions = [];
    successful_generations = [];
    
    % --- 1. FILTRAGEM DOS PACOTES (LÓGICA DE PREEMPÇÃO) ---
    % Um pacote k só termina se: (Chegada_k + Serviço_k) <= Chegada_k+1
    % Caso contrário, ele é interrompido pela chegada de k+1.
    
    for k = 1:(N-1)
        potential_finish_time = generation_times(k) + service_times(k);
        next_arrival_time = generation_times(k+1);
        
        if potential_finish_time <= next_arrival_time
            % O pacote conseguiu terminar antes do próximo chegar
            successful_completions = [successful_completions; potential_finish_time];
            successful_generations = [successful_generations; generation_times(k)];
        else
            % O pacote foi preemptado (morto). Não contribui para o AoI.
            % (Não fazemos nada, apenas ignoramos ele no cálculo da área)
        end
    end
    
    % O último pacote sempre termina (pois não há k+1 para interrompê-lo)
    % Assumindo que a simulação permite que o último termine.
    last_finish = generation_times(N) + service_times(N);
    successful_completions = [successful_completions; last_finish];
    successful_generations = [successful_generations; generation_times(N)];
    
    % --- 2. CÁLCULO DA ÁREA SOB A CURVA ---
    area_total = 0;
    completion_time_prev = 0;
    generation_time_prev = 0; % Assume-se sistema vazio em t=0
    
    num_success = length(successful_completions);
    
    for k = 1:num_success
        C_curr = successful_completions(k);
        G_curr = successful_generations(k);
        
        % Intervalo de tempo entre a atualização anterior e a atual
        interval = C_curr - completion_time_prev;
        
        % Altura do triângulo/trapézio
        % H1: Altura do AoI logo após a atualização anterior (t_prev - g_prev)
        H1 = completion_time_prev - generation_time_prev;
        
        % H2: Altura do AoI imediatamente antes da nova atualização (t_curr - g_prev)
        % Nota: O AoI cresce linearmente com inclinação 1
        H2 = C_curr - generation_time_prev;
        
        % Área do trapézio
        area_step = 0.5 * (H1 + H2) * interval;
        area_total = area_total + area_step;
        
        % Atualiza para o próximo passo
        completion_time_prev = C_curr;
        generation_time_prev = G_curr; % O AoI cai instantaneamente para (C_curr - G_curr)
    end
    
    % Média do AoI
    % Se nenhum pacote passou (caso extremo), retorna 0 ou NaN para evitar divisão por zero
    if completion_time_prev == 0
        avg_aoi = 0; 
    else
        avg_aoi = area_total / completion_time_prev;
    end
end