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

