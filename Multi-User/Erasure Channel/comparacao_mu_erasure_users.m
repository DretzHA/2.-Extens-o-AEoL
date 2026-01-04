function comparacao_mu_erasure_users()
    clc; clear; close all;

    %% CONFIGURAÇÕES GERAIS
    sim_params.k = 50;               % Tamanho do pacote (símbolos)
    sim_params.delta = 0.5;          % Probabilidade FIXA de apagamento
    sim_params.T_max = 200000;       % Tempo de simulação (slots)
    sim_params.lambda_per_user = 1/200; % Taxa de chegada POR USUÁRIO
    
    % Variação do número de usuários para o eixo X
    user_counts = [1, 3, 5, 10, 15, 20, 25, 30]; 
    
    % Armazenar resultados
    aoi_preemptive = zeros(size(user_counts));
    aoi_non_preemptive = zeros(size(user_counts));

    fprintf('--- Simulação Multi-User com Erasure Channel ---\n');
    fprintf('Delta (Prob. Erro) fixo em: %.2f\n', sim_params.delta);
    fprintf('Taxa por usuário (lambda): %.4f\n\n', sim_params.lambda_per_user);
    fprintf('%-10s %-20s %-20s\n', 'N Users', 'AoI (Preemptive)', 'AoI (Non-Preempt)');

    %% LOOP PRINCIPAL (Varia N Usuários)
    for i = 1:length(user_counts)
        N = user_counts(i);
        sim_params.N_users = N;
        
        % O lambda total do sistema aumenta com o número de usuários
        % Mas a simulação lida com gerações independentes ou um stream agregado.
        % Vamos passar o N para dentro da função.
        
        % 1. Simulação Com Preempção (Global LCLS)
        aoi_preemptive(i) = run_multiuser_simulation(sim_params, 'Preemptive');
        
        % 2. Simulação Sem Preempção (Waiting / Buffer=1)
        aoi_non_preemptive(i) = run_multiuser_simulation(sim_params, 'Non-Preemptive');
        
        fprintf('%-10d %-20.4f %-20.4f\n', N, aoi_preemptive(i), aoi_non_preemptive(i));
    end

    %% PLOTAGEM DOS RESULTADOS
    plot_results(user_counts, aoi_preemptive, aoi_non_preemptive, sim_params);
end

%% FUNÇÃO DE SIMULAÇÃO MULTI-USER
function mean_avg_aoi = run_multiuser_simulation(params, policy)
    % Parâmetros
    N = params.N_users;
    lambda = params.lambda_per_user;
    k = params.k;
    delta = params.delta;
    T_max = params.T_max;
    
    % Estado dos Usuários
    % Cada usuário tem seu próprio Age atual e tempo da última atualização
    ages = zeros(1, N);       % Age instantâneo de cada usuário (começa em 0 idealmente)
    last_generation_times = zeros(1, N); % Data de nascimento do pacote entregue
    
    % Para calcular a área (AoI médio) de cada usuário
    area_aoi_users = zeros(1, N);
    
    % Estado do Servidor
    current_time = 0;
    is_busy = false;
    service_finish_time = inf;
    
    % Pacote sendo servido atualmente: [user_id, generation_time]
    packet_in_service = []; 
    
    % Fila de Espera (Apenas para Non-Preemptive com Buffer=1)
    % Guardamos apenas UM pacote na espera (o mais fresco globalmente ou por usuário?)
    % Modelo comum em MU-AoI: Buffer compartilhado de tamanho 1 (substituição)
    queue_packet = []; % [user_id, generation_time]
    
    % Geração inicial de chegadas: uma para cada usuário
    % Event list: [time, type, user_id]
    % type: 1 = Arrival, 2 = Departure
    next_arrivals = zeros(1, N);
    for u = 1:N
        next_arrivals(u) = exprnd(1/lambda);
    end
    
    %% LOOP DE EVENTOS DISCRETOS
    while current_time < T_max
        % 1. Encontrar o próximo evento (Chegada de alguém ou Fim de serviço)
        [min_arr_time, u_idx] = min(next_arrivals);
        
        if min_arr_time < service_finish_time
            event_time = min_arr_time;
            event_type = 'ARRIVAL';
            event_user = u_idx;
        else
            event_time = service_finish_time;
            event_type = 'DEPARTURE';
            event_user = 0; % Definido pelo pacote em serviço
        end
        
        dt = event_time - current_time;
        current_time = event_time;
        
        % 2. Atualizar Áreas de AoI para TODOS os usuários
        % O Age de cada usuário cresce linearmente com o tempo
        for u = 1:N
            % Área trapézio: (Age_old + Age_new)*dt/2
            % Age_new = Age_old + dt
            current_age = current_time - last_generation_times(u);
            prev_age = current_age - dt;
            
            area_aoi_users(u) = area_aoi_users(u) + (prev_age + current_age) * dt / 2;
        end
        
        % 3. Tratar Evento
        switch event_type
            case 'ARRIVAL'
                % Nova chegada para o usuário 'event_user'
                gen_time = current_time;
                
                if ~is_busy
                    % Servidor livre: atende imediatamente
                    is_busy = true;
                    packet_in_service = [event_user, gen_time];
                    
                    % Gera tempo de serviço estocástico (Erasure)
                    s_time = generate_service_time(k, delta);
                    service_finish_time = current_time + s_time;
                    
                else
                    % Servidor Ocupado
                    if strcmp(policy, 'Preemptive')
                        % PREEMPÇÃO GLOBAL (LCLS)
                        % O pacote novo (de qualquer usuário) chuta o atual
                        packet_in_service = [event_user, gen_time];
                        
                        % Resetar tempo de serviço (novo canal, nova tentativa)
                        s_time = generate_service_time(k, delta);
                        service_finish_time = current_time + s_time;
                        
                    else
                        % NON-PREEMPTIVE (Com Buffer de Troca)
                        % O pacote vai para a espera. Se já tinha um na espera, substitui.
                        % (Buffer Size = 1 global, política comum para maximizar frescor em fila)
                        queue_packet = [event_user, gen_time];
                    end
                end
                
                % Agendar próxima chegada para este usuário específico
                next_arrivals(event_user) = current_time + exprnd(1/lambda);
                
            case 'DEPARTURE'
                % Pacote terminou de ser enviado
                finished_user = packet_in_service(1);
                finished_gen_time = packet_in_service(2);
                
                % Atualiza a "última geração recebida" desse usuário
                % O Age instantâneo cai para (CurrentTime - GenTime)
                % Nota: Só atualizamos se esse pacote for mais novo do que o que o usuário já tem
                if finished_gen_time > last_generation_times(finished_user)
                    last_generation_times(finished_user) = finished_gen_time;
                end
                
                is_busy = false;
                packet_in_service = [];
                service_finish_time = inf;
                
                % Se há pacote na fila (apenas Non-Preemptive usa isso)
                if ~isempty(queue_packet)
                    is_busy = true;
                    packet_in_service = queue_packet;
                    queue_packet = []; % Esvazia fila
                    
                    s_time = generate_service_time(k, delta);
                    service_finish_time = current_time + s_time;
                end
        end
    end
    
    % Retorna a média das médias dos usuários
    mean_aoi_per_user = area_aoi_users / T_max;
    mean_avg_aoi = mean(mean_aoi_per_user);
end

%% GERADOR DE TEMPO DE SERVIÇO (Binomial Negativa)
function S = generate_service_time(k, delta)
    % Retorna o número total de slots (sucessos + falhas)
    % para obter k sucessos com probabilidade (1-delta)
    prob_sucesso = 1 - delta;
    
    if prob_sucesso >= 1 || delta == 0
        S = k;
    else
        % nbinrnd retorna o número de FALHAS antes de k sucessos
        n_failures = nbinrnd(k, prob_sucesso);
        S = k + n_failures; 
    end
end

function plot_results(users, aoi_p, aoi_np, params)
    figure;
    plot(users, aoi_p, '-o', 'LineWidth', 2, 'DisplayName', 'Preemptive (LCLS)');
    hold on;
    plot(users, aoi_np, '-s', 'LineWidth', 2, 'DisplayName', 'Non-Preemptive (Wait)');
    
    xlabel('Número de Usuários (N)');
    ylabel('Average AoI por Usuário (Slots)');
    title(sprintf('Comparação AoI Multi-User com Erasure (\\delta=%.2f, k=%d)', params.delta, params.k));
    legend('Location', 'NorthWest');
    grid on;
    
    % Adicionar anotação de parâmetros
    dim = [.65 .5 .3 .3];
    str = {sprintf('Lambda/User: %.4f', params.lambda_per_user), ...
           sprintf('Prob. Apagamento: %.2f', params.delta)};
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');
end