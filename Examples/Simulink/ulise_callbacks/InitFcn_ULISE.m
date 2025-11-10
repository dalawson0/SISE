% File: InitFcn_ULISE.m
%
% Description:
%   Initialization script for the ULISE Simulink model.
%   Loads state-space matrices and Parameter into the model workspace.
%   Runs automatically via the model’s InitFcn callback.
%
% Author: Daniel Lawson
% Date:   October 10, 2025

%% Simulation Setup Parameter

SimSetup = struct();

SimSetup.Ts = 1e0; % Sampling time
SimSetup.K= 1e3; % Total number of time steps for simulation
SimSetup.T = SimSetup.Ts * SimSetup.K; % total simulation time
SimSetup.t  = (0:SimSetup.K-1)' * SimSetup.Ts; % Time vector
%% Define Plant Parameter


A = [ 0.5    2   0   0   0;
    0      0.2 1   0   1; 
    0      0   0.3 0   1; 
    0      0   0   0.7 1;
    0      0   0   0   0.1 ]; 

% B = zeros(size(A,1),1);
% 
% C = blkdiag(1,1,1,1,1);
% 
% D = zeros(size(C,1),1);
% 
% G = [ 1  0   -0.3;
%     1   0   0;
%     0   0   0;
%     0   0   0;
%     0   0   0] ;
% p = size(G, 2);
% 
% H = [ 0  0   1;
%     0  0   0; 
%     0  1   0; 
%     0  0   0; 
%     0  0   0 ]; 

% x0 = zeros(size(A,1),1);
% 
% R = 1e-2 * [  1    0   0   0.5 0;
%             0    1   0   0   0.3; 
%             0    0   1   0   0; 
%             0.5  0   0   1   0; 
%             0    0.3 0   0   1  ];
% 
% Q = 1e-4 * [  1    0   0   0   0; 
%             0    1   0.5 0   0; 
%             0    0.5 1   0   0; 
%             0    0   0   1   0; 
%             0    0   0   0   1  ]; 
% 
% 
xhat0 = zeros(size(A,1),1);
P_x0 = zeros(size(A));

% % Transmission zeros
% t_zeros = tzero(A,G,C,H);
% if rank([A G; C H]) ~= size(A,1)+size(G,2)
%     error('Error. System is not strongly detectable.')
% end
% 
% for i = 1:length(t_zeros)
%     if abs(t_zeros(i)) > 1
%         error('Error. System is not strongly detectable')
% 
%     end
% end

% Number of unkown inputs
p = 3;

% Numbe of known inputs
m = 0;

% Number of states
n = 5;

% Number of Outputs
l = n;

%% Generate the process noise and measurement noise vectors
% 
Q = 1e-4 * [  1    0   0   0   0; 
            0    1   0.5 0   0; 
            0    0.5 1   0   0; 
            0    0   0   1   0; 
            0    0   0   0   1  ]; 

R = 1e-2 * [  1    0   0   0.5 0;
            0    1   0   0   0.3; 
            0    0   1   0   0; 
            0.5  0   0   1   0; 
            0    0.3 0   0   1  ];

SigmaQ = chol(Q);
w = (randn(SimSetup.K,n)*SigmaQ)';

SigmaR = chol(R);
v = (randn(SimSetup.K,l)*SigmaR)';

%% Simulate 'true' known input, u, and unknown input, d. 

u = zeros(m,SimSetup.K);

d = zeros(p,SimSetup.K); 
d(1,500:700) = 1*ones(size(d(1,500:700)));

for j = 100:800
    d(2,j) = 1/700*(j-100);
end

% Modify the third disturbance component
d(3,[500:549 600:649 700:749]) = 3*ones(size(d(3,[500:549 600:649 700:749])));
d(3,[550:599 650:699 750:799]) = -3*ones(size(d(3,[550:599 650:699 750:799])));

%% Make timeseries input data 

InputData = struct();

% Add known and unknown input data to struct
InputData.u=setinterpmethod(timeseries(u, SimSetup.t),'zoh');
InputData.u.Name = "timeseries known input";

InputData.d=setinterpmethod(timeseries(d, SimSetup.t),'zoh');
InputData.d.Name = "timeseries unknown input";

% Add process and measurement noise timeseries to struct
InputData.w = setinterpmethod(timeseries(w, SimSetup.t),'zoh'); 
InputData.w.Name = "timeseries process noise";

InputData.v = setinterpmethod(timeseries(v, SimSetup.t),'zoh');
InputData.v.Name = "timeseries measurement noise";
