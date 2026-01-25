function [avg_aoi_per_user, avg_aoi_system] = queue_logic_multi(arrivals, services, num_users, type)

    if strcmp(type, 'NP')
        completed = run_lcfs_np(arrivals, services);
    else
        completed = run_lcfs_p(arrivals, services);
    end
    
    % =====================================================================
    % 1. System AoI
    if isempty(completed)
        avg_aoi_system = 0;
    else
        data_sys = completed(:, 1:2); % [Completion, Generation]
        avg_aoi_system = calc_aoi_area(data_sys);
    end

    % =====================================================================
    % 2. User AoI
    % =====================================================================
    avg_aoi_per_user = zeros(num_users, 1);
    
    for u = 1:num_users
        mask = (completed(:, 3) == u);
        data_u = completed(mask, 1:2);
        
        if isempty(data_u)
            avg_aoi_per_user(u) = 0;
        else
            avg_aoi_per_user(u) = calc_aoi_area(data_u);
        end
    end
end


function avg_val = calc_aoi_area(data)
    % data: [Completion_Time, Generation_Time]
    
    % Sort by time
    data = sortrows(data, 1);
    C = data(:, 1); 
    G = data(:, 2); 
    
    area = 0; 
    t_now = 0;      
    current_aoi_start = 0; 
    
    for k = 1:length(C)
        t_fin = C(k);
        t_gen = G(k);
        
        if t_gen < current_aoi_start

        end

        % Interval [t_now, t_fin]
        dt = t_fin - t_now;
        
        % Area
        age_start = t_now - current_aoi_start;
        age_end   = t_fin - current_aoi_start;
        
        area = area + (age_start + age_end) * dt / 2;
        
        t_now = t_fin;      
        current_aoi_start = t_gen; 
    end
    
    if t_now > 0
        avg_val = area / t_now;
    else
        avg_val = 0;
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
            
        end
    end
    
    last_idx = N;
    finish_last = arr(last_idx, 1) + srv(last_idx);
    completed = [completed; finish_last, arr(last_idx, 1), arr(last_idx, 2)];
end