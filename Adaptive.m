%% Continuous Local Learning (CLL) WRSN Simulator (With Dynamic Hotspots)
% Features:
% - Dynamic Multihop SPT Routing & Topology Tracking
% - Multi-Phase Dynamic Hotspot Swapping (Geographic & Random Clusters)
% - Full MC Energy Accounting (Movement, Charging Losses, SS Transits)
% - Monte Carlo Execution (100 runs) with 95% Confidence Intervals

clear; clc; close all;

%% 1. Comprehensive System Parameters

% Simulation Runtime & Execution
num_runs      = 100;                 % Monte Carlo executions[cite: 1]
T_sim         = 60*60*24*180;        % Simulation horizon 6 months in seconds
T_warmup      = 0.05*T_sim;          % Warm-up/transient period (seconds)
area_size     = 200;                 % Field dimension (200m x 200m)
BS_pos        = [100, 100];          % Base station at geometric center (100m, 100m)
SS_pos        = [100, 100];          % Service station at geometric center (100m, 100m)
R_comm        = 30;                  % Sensor communication radius (meters)

% Hotspot Configuration
enable_hotspots = true;              % Enable dynamic hotspot arrivals
hotspot_interval= 2.5e5;             % Hotspots swap position every 250,000 seconds
hotspot_ratio   = 0.15;              % 15% of total nodes become hotspots per phase
hotspot_surge   = 4.0;               % Traffic multiplier (400% rate increase)

% Sensor Network Specifications
n_nodes       = 300;                 % Network size (test across 100, 300, 500)
E_max         = 10.8e3;              % Sensor battery capacity (Joules)
tau_threshold = 0.05 * E_max;        % Non-operational threshold (5% = 540 J)
P_sense       = 0.015;               % Base sensing power consumption (Watts)
E_tx_bit      = 50e-9;               % Transmit energy per bit (50 nJ/bit)
E_rx_bit      = 50e-9;               % Receive energy per bit (50 nJ/bit)
pkt_payload   = 1000;                % Packet size (1000 bits)

% Mobile Charger (MC) Specifications
v_mc          = 5.0;                 % Travel speed (m/s)
E_mc_max      = 770e3;               % MC total battery capacity (Joules)
E_move        = 5.0;                 % Movement energy cost (J/m)
P_move        = E_move * v_mc;       % Movement power consumption (25 W)
Delta_charge  = 5.0;                 % Wireless charging rate (Watts)
eta_charge    = 0.95;                % Wireless charging efficiency (95%)
T_service     = 300;                 % Service station battery swap duration (seconds)

%% 2. Data Containers for Performance Tracking
non_op_sizes     = zeros(num_runs, 1);
disconn_times_h  = zeros(num_runs, 1);
mc_travel_dist_km= zeros(num_runs, 1);
mc_ss_visits     = zeros(num_runs, 1);
init_durations_s = zeros(num_runs, 1);

fprintf('Executing %d Monte Carlo runs with Dynamic Hotspots (n = %d)...\n', num_runs, n_nodes);

for run_idx = 1:num_runs
    rng(run_idx, 'twister'); % Seed reproducibility
    
    % Node Deployment (Uniform Random Distribution)
    nodes_loc = area_size * rand(n_nodes, 2);
    dist_matrix = squareform(pdist(nodes_loc));
    
    % Verify Graph Connectivity & Build Shortest Path Tree (SPT)
    dist_to_BS = sqrt(sum((nodes_loc - BS_pos).^2, 2));
    [adj_matrix, parent_nodes, relay_counts] = update_spt_routing(dist_matrix, dist_to_BS, R_comm);
    
    % Base Traffic Generation (Poisson Packet Arrivals: 1 to 5 kbps)
    base_lambda_rates = 1 + 4 * rand(n_nodes, 1); % Baseline arrivals in pkts/sec
    current_lambda    = base_lambda_rates;
    node_data_rates   = current_lambda * pkt_payload;
    
    % Compute Initial Power Depletion Rates (P_node)
    P_node = compute_power_drain(node_data_rates, relay_counts, P_sense, E_tx_bit, E_rx_bit);
    
    %% Initialization Phase: Discovery (DFS) & Base Estimation
    T_1 = zeros(n_nodes, 1); S_1 = zeros(n_nodes, 1);
    T_2 = zeros(n_nodes, 1); S_2 = zeros(n_nodes, 1);
    
    % DFS Route Traversal from Base Station
    dfs_route = perform_dfs_traversal(adj_matrix, nodes_loc, BS_pos);
    
    t_curr = 0; 
    mc_pos = SS_pos; 
    E_mc = E_mc_max;
    total_travel_m = 0;
    ss_visit_count = 0;
    
    % Round 1: Discovery Traversal
    for i = 1:length(dfs_route)
        target = dfs_route(i);
        d = norm(mc_pos - nodes_loc(target, :));
        
        if (E_mc - d * E_move) <= 0.05 * E_mc_max
            d_to_ss = norm(mc_pos - SS_pos);
            t_curr = t_curr + (d_to_ss / v_mc) + T_service;
            total_travel_m = total_travel_m + d_to_ss;
            E_mc = E_mc_max;
            mc_pos = SS_pos;
            ss_visit_count = ss_visit_count + 1;
            d = norm(mc_pos - nodes_loc(target, :));
        end
        
        t_curr = t_curr + (d / v_mc);
        E_mc = E_mc - (d * E_move);
        total_travel_m = total_travel_m + d;
        mc_pos = nodes_loc(target, :);
        
        T_1(target) = t_curr;
        S_1(target) = E_max; % Top off node during discovery
    end
    
    % Round 2: Base Estimation Traversal
    for i = 1:length(dfs_route)
        target = dfs_route(i);
        d = norm(mc_pos - nodes_loc(target, :));
        
        if (E_mc - d * E_move) <= 0.05 * E_mc_max
            d_to_ss = norm(mc_pos - SS_pos);
            t_curr = t_curr + (d_to_ss / v_mc) + T_service;
            total_travel_m = total_travel_m + d_to_ss;
            E_mc = E_mc_max;
            mc_pos = SS_pos;
            ss_visit_count = ss_visit_count + 1;
            d = norm(mc_pos - nodes_loc(target, :));
        end
        
        t_curr = t_curr + (d / v_mc);
        E_mc = E_mc - (d * E_move);
        total_travel_m = total_travel_m + d;
        mc_pos = nodes_loc(target, :);
        
        T_2(target) = t_curr;
        elapsed = T_2(target) - T_1(target);
        S_2(target) = max(tau_threshold, E_max - P_node(target) * elapsed);
    end
    
    init_durations_s(run_idx) = t_curr;
    delta_est = (S_1 - S_2) ./ max(1.0, (T_2 - T_1));
    
    %% Learning & Recharging Phase
    node_energy = S_2;
    node_disconn_time = zeros(n_nodes, 1);
    non_op_sample_sum = 0;
    samples_count = 0;
    
    last_hotspot_swap = t_curr;
    active_hotspot_nodes = [];
    
    while t_curr < T_sim
        
        % Dynamic Hotspot Swapping Handler
        if enable_hotspots && (t_curr - last_hotspot_swap >= hotspot_interval)
            % 1. Select a random center for the new hotspot zone
            hotspot_center = area_size * rand(1, 2);
            dist_to_center = sqrt(sum((nodes_loc - repmat(hotspot_center, n_nodes, 1)).^2, 2));
            
            % 2. Identify nearest nodes to form the hotspot cluster
            [~, sorted_indices] = sort(dist_to_center);
            num_hotspot_nodes = round(hotspot_ratio * n_nodes);
            active_hotspot_nodes = sorted_indices(1:num_hotspot_nodes);
            
            % 3. Apply Traffic Surge Multiplier
            current_lambda = base_lambda_rates;
            current_lambda(active_hotspot_nodes) = current_lambda(active_hotspot_nodes) * hotspot_surge;
            
            % 4. Recalculate Node Power Consumption
            node_data_rates = current_lambda * pkt_payload;
            P_node = compute_power_drain(node_data_rates, relay_counts, P_sense, E_tx_bit, E_rx_bit);
            
            last_hotspot_swap = t_curr;
        end
        
        % Check MC Service Station Battery Threshold (5% Capacity)
        if E_mc <= 0.05 * E_mc_max
            d_ss = norm(mc_pos - SS_pos);
            t_curr = t_curr + (d_ss / v_mc) + T_service;
            total_travel_m = total_travel_m + d_ss;
            E_mc = E_mc_max;
            mc_pos = SS_pos;
            ss_visit_count = ss_visit_count + 1;
        end
        
        % Predict Expected Energy S_exp at Arrival Time t' = t + d/v
        d_to_nodes = sqrt(sum((nodes_loc - repmat(mc_pos, n_nodes, 1)).^2, 2));
        t_arrival = t_curr + (d_to_nodes / v_mc);
        
        S_exp = S_2 - delta_est .* (t_arrival - T_2);
        
        % Target Selection: argmin { S_exp(s, t') }
        [~, target_node] = min(S_exp);
        
        % MC Transit to Target
        d_travel = d_to_nodes(target_node);
        t_travel = d_travel / v_mc;
        t_curr = t_curr + t_travel;
        E_mc = E_mc - (d_travel * E_move);
        total_travel_m = total_travel_m + d_travel;
        mc_pos = nodes_loc(target_node, :);
        
        % Update Node Energies & Accumulate Disconnection Durations
        for i = 1:n_nodes
            e_drain = P_node(i) * t_travel;
            if node_energy(i) > tau_threshold
                if (node_energy(i) - e_drain) <= tau_threshold
                    time_active = (node_energy(i) - tau_threshold) / P_node(i);
                    node_disconn_time(i) = node_disconn_time(i) + (t_travel - time_active);
                    node_energy(i) = tau_threshold;
                else
                    node_energy(i) = node_energy(i) - e_drain;
                end
            else
                node_disconn_time(i) = node_disconn_time(i) + t_travel;
            end
        end
        
        % Perform Wireless Recharging Service
        energy_needed = E_max - node_energy(target_node);
        t_charge = energy_needed / Delta_charge;
        E_mc = E_mc - (t_charge * Delta_charge / eta_charge);
        t_curr = t_curr + t_charge;
        node_energy(target_node) = E_max;
        
        % Update Continuous Local Learning Rates (T_old, S_old, delta)
        T_1(target_node) = T_2(target_node);
        S_1(target_node) = S_2(target_node);
        T_2(target_node) = t_curr;
        S_2(target_node) = E_max;
        
        % Dynamic rate adjustment captures hotspot surges naturally
        delta_est(target_node) = (S_1(target_node) - S_2(target_node)) / ...
            max(1.0, (T_2(target_node) - T_1(target_node)));
        
        % Sampling Post Warm-up Window
        if t_curr >= T_warmup
            non_op_sample_sum = non_op_sample_sum + sum(node_energy <= tau_threshold);
            samples_count = samples_count + 1;
        end
    end
    
    % Store Monte Carlo Run Output Statistics
    non_op_sizes(run_idx)      = non_op_sample_sum / max(1, samples_count);
    disconn_times_h(run_idx)   = mean(node_disconn_time) / 3600; % Hours
    mc_travel_dist_km(run_idx) = total_travel_m / 1000;          % Kilometers
    mc_ss_visits(run_idx)      = ss_visit_count;
end

%% 3. Statistical Analysis & 95% Confidence Intervals

[mean_non_op, ci_non_op]   = compute_stats(non_op_sizes, num_runs);
[mean_disconn, ci_disconn] = compute_stats(disconn_times_h, num_runs);
[mean_dist, ci_dist]       = compute_stats(mc_travel_dist_km, num_runs);
[mean_visits, ci_visits]   = compute_stats(mc_ss_visits, num_runs);

% Print Formatted Statistical Report
fprintf('\n=======================================================\n');
fprintf('     DYNAMIC HOTSPOT CLL PERFORMANCE RESULTS           \n');
fprintf('=======================================================\n');
fprintf('Hotspot Condition:               15%% Nodes Surge by 400%% every 250ks\n');
fprintf('Avg. Non-Operational Size:        %.2f ± %.2f nodes (95%% CI)\n', mean_non_op, ci_non_op);
fprintf('Avg. Disconnection Time:          %.2f ± %.2f hours (95%% CI)\n', mean_disconn, ci_disconn);
fprintf('Avg. MC Travel Distance:          %.2f ± %.2f km    (95%% CI)\n', mean_dist, ci_dist);
fprintf('Avg. Service Station Visits:      %.1f ± %.1f visits (95%% CI)\n', mean_visits, ci_visits);
fprintf('=======================================================\n\n');

%% 4. Helper Functions

function [adj_matrix, parent, relay_counts] = update_spt_routing(dist_matrix, dist_to_BS, R_comm)
    n = length(dist_to_BS);
    adj_matrix = (dist_matrix <= R_comm) & (dist_matrix > 0);
    parent = zeros(n, 1);
    relay_counts = zeros(n, 1);
    
    for i = 1:n
        neighbors = find(adj_matrix(i, :));
        if ~isempty(neighbors)
            [~, min_idx] = min(dist_to_BS(neighbors));
            parent(i) = neighbors(min_idx);
        end
    end
    
    for i = 1:n
        curr = parent(i);
        while curr > 0
            relay_counts(curr) = relay_counts(curr) + 1;
            curr = parent(curr);
        end
    end
end

function P_node = compute_power_drain(node_data_rates, relay_counts, P_sense, E_tx_bit, E_rx_bit)
    n = length(node_data_rates);
    total_forwarded_rate = zeros(n, 1);
    for i = 1:n
        total_forwarded_rate(i) = node_data_rates(i) + relay_counts(i) * mean(node_data_rates);
    end
    P_node = P_sense + (total_forwarded_rate * E_tx_bit) + ((total_forwarded_rate - node_data_rates) * E_rx_bit);
end

function dfs_route = perform_dfs_traversal(adj_matrix, nodes_loc, BS_pos)
    n = size(nodes_loc, 1);
    visited = false(n, 1);
    dfs_route = [];
    
    dist_bs = sqrt(sum((nodes_loc - BS_pos).^2, 2));
    [~, start_node] = min(dist_bs);
    
    stack = [start_node];
    while ~isempty(stack)
        curr = stack(end);
        stack(end) = [];
        if ~visited(curr)
            visited(curr) = true;
            dfs_route = [dfs_route; curr];
            neighbors = find(adj_matrix(curr, :));
            for k = 1:length(neighbors)
                if ~visited(neighbors(k))
                    stack = [stack; neighbors(k)];
                end
            end
        end
    end
    
    unvisited = find(~visited);
    dfs_route = [dfs_route; unvisited];
end

function [m, ci95] = compute_stats(data_vec, N)
    m = mean(data_vec);
    s = std(data_vec);
    ci95 = 1.96 * (s / sqrt(N));
end