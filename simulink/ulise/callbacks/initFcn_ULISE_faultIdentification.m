% File: InitFcn_ULISE.m
%
% Description:
%   Initialization script for the ULISE Simulink model. Defines simulation
%   parameters and loads system matrices into the model workspace. It also
%   generates custom known and unknown inputs and stroes them in a struct
%   for the simulation. Runs automatically via the model’s InitFcn
%   callback. 
%
%   Authorship:
%       D. Lawson, S. Macero, E. Gah, and S.Z. Yong (2025). SISE Toolbox
%       for MATLAB (https://github.com/dalawson0/SISE.git), GitHub. 

%% Simulation Setup Parameters

SimSetup = struct();

SimSetup.Ts = 1e0; % Sampling time
SimSetup.K  = 1e3; % Total number of time steps for simulation

%% Define Plant Parameters

System = struct();

% State Trasition Matrix
System.A = [0.5  2    0    0    0;
            0    0.2  1    0    1; 
            0    0    0.3  0    1; 
            0    0    0    0.7  1;
            0    0    0    0    0.1]; 

% Control Input Matrix
System.B = zeros(5, 0);

% Measurement Matrix
System.C = blkdiag(1,1,1,1,1);

% Measurement feedthrough from u
System.D = zeros(5,0);

% Disturbance input matrix
System.G = [1   0   -0.3;
            1   0   0;
            0   0   0;
            0   0   0;
            0   0   0];

% Disturbance feedthrough to output
System.H = [0  0  1;
            0  0  0; 
            0  1  0; 
            0  0  0; 
            0  0  0]; 

% Process noise covariance
System.Q = 1e-4 * [1    0    0    0   0; 
                   0    1    0.5  0   0; 
                   0    0.5  1    0   0; 
                   0    0    0    1   0; 
                   0    0    0    0   1]; 


% Measurement noise covariance
System.R = 1e-2 * [1    0    0   0.5  0;
                   0    1    0   0    0.3; 
                   0    0    1   0    0; 
                   0.5  0    0   1    0; 
                   0    0.3  0   0    1];

% Initial state vector
System.x0 = zeros(5,1);

% Initial estimate of the state vector
System.xhat0 = zeros(5,1);

% Initial state covariance matrix
System.P_x0 = eye(size(System.A,1))*1e6;

%% Simulate 'true' known input, u, and unknown input, d (change if desired)

u = zeros(0, SimSetup.K);

d = zeros(3, SimSetup.K);
d(1,500:700) = 1*ones(size(d(1,500:700)));

for j = 100:800
    d(2,j) = 1/700*(j-100);
end

% Modify the third disturbance component
d(3,[500:549 600:649 700:749]) = 3*ones(size(d(3,[500:549 600:649 700:749])));
d(3,[550:599 650:699 750:799]) = -3*ones(size(d(3,[550:599 650:699 750:799])));

%% Make timeseries input data and update known input matrices 

% Default filename (change if desired, but remember to also change the name in the Signal Editor block)
System.filename = 'ULISE_default_signals.mat'; 

% Internal processing (do not change)
[System.B, System.D, input_data] = preprocessing_ULISE(System.A, System.B, ...
    System.C, System.D, System.G, SimSetup.K, SimSetup.Ts, u, d, System.filename);

