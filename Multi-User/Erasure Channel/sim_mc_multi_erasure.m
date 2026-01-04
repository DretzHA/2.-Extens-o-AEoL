function [Sim_NP, Sim_P] = sim_mc_multi_erasure(conf)

    num_L = length(conf.lambda_total_vec);
    Sim_NP = zeros(conf.num_users, num_L);
    Sim_P  = zeros(conf.num_users, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        
        % Follwing the Eq. 1 - The Age of Information: Real-Time Status Updating by Multiple Sources
        % we divide the total lambda equally by the number of users,
        % maintening the relation of rho = sum(rho_i)

        lambda_user = lambda_tot / conf.num_users;
        
        acc_aoi_np = zeros(conf.num_users, 1);
        acc_aoi_p  = zeros(conf.num_users, 1);
        
        for k_mc = 1:conf.N_MC
            % Gerenerate the arrival time of the updates, following poisson
            % process
            all_arrivals = [];
            for u = 1:conf.num_users
                dt = exprnd(1/lambda_user, conf.N_updates, 1); 
                t_arr = cumsum(dt); 
                all_arrivals = [all_arrivals; t_arr, ones(conf.N_updates, 1)*u];
            end
            
            all_arrivals = sortrows(all_arrivals, 1);
            
            % Services Time (Erasure Channel)
            n_pkts = size(all_arrivals, 1);
            services_needed = get_erasure_service(n_pkts, conf.k, conf.delta);
            services = services_needed* conf.Symbol_Duration;;
            
            % Process Queues
            
            % LCFS Non-Preemptive (Preemptive in Waiting)
            aoi_np = queue_logic_multi(all_arrivals, services, conf.num_users, 'NP');
            acc_aoi_np = acc_aoi_np + aoi_np;
            
            % LCFS Preemptive (Preemptice in Service)
            aoi_p = queue_logic_multi(all_arrivals, services, conf.num_users, 'P');
            acc_aoi_p = acc_aoi_p + aoi_p;
        end
        
        % Mean Values
        mean_aoi_np = acc_aoi_np / conf.N_MC;
        mean_aoi_p  = acc_aoi_p  / conf.N_MC;
        
        for u = 1:conf.num_users
            v = conf.velocities(u);
            Sim_NP(u, j) = v * mean_aoi_np(u);
            Sim_P(u, j)  = v * mean_aoi_p(u);
        end
    end
end

%% FUNCTION ERASURE CHANNEL
function S = get_erasure_service(n_packets, k, delta)

    prob_sucess = 1 - delta;
    
    if prob_sucess >= 1 || delta == 0
        S = ones(n_packets, 1) * k;
    else
        % nbinrnd(k, p) returns the number of fails before k sucess
        n_failures = nbinrnd(k, prob_sucess, n_packets, 1);
        
        % Total Time (per time unit) = k + Fails
        S = k + n_failures;
    end
end