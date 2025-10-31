% File: InitFcn_ULISE.m
%
% Description:
%   Initialization script for the ULISE Simulink model.
%   Loads state-space matrices and Parameter into the model workspace.
%   Runs automatically via the model’s InitFcn callback.
%
% Author: Daniel Lawson
% Date:   October 10, 2025

%% Define Plant Parameter


A = [ 0.5    2   0   0   0;
    0      0.2 1   0   1; 
    0      0   0.3 0   1; 
    0      0   0   0.7 1;
    0      0   0   0   0.1 ]; 

B = zeros(size(A,1),1);

C = blkdiag(1,1,1,1,1);

D = zeros(size(C,1),1);

G = [ 1  0   -0.3;
    1   0   0;
    0   0   0;
    0   0   0;
    0   0   0] ;
p = size(G, 2);

H = [ 0  0   1;
    0  0   0; 
    0  1   0; 
    0  0   0; 
    0  0   0 ]; 

x0 = zeros(size(A,1),1);

R = 1e-2 * [  1    0   0   0.5 0;
            0    1   0   0   0.3; 
            0    0   1   0   0; 
            0.5  0   0   1   0; 
            0    0.3 0   0   1  ];

Q = 1e-4 * [  1    0   0   0   0; 
            0    1   0.5 0   0; 
            0    0.5 1   0   0; 
            0    0   0   1   0; 
            0    0   0   0   1  ]; 

% Singular Value Decomposition of H
r = rank(H);
[U, S, V] = svd(H);

U1 = U(:, 1:r); 
U2 = U(:, r+1:end);

Sigma = S(1:r, 1:r);
 
V1 = V(:, 1:r);
V2 = V(:, r+1:end);

xhat0 = zeros(size(A,1),1);

% T1 = U1' - U1'*R*U2*inv(U2'*R*U2)*U2';
% T2 = U2';

% Transmission zeros
t_zeros = tzero(A,G,C,H);
if rank([A G; C H]) ~= size(A,1)+size(G,2)
    error('Error. System is not strongly detectable.')
end

for i = 1:length(t_zeros)
    if abs(t_zeros(i)) > 1
        error('Error. System is not strongly detectable')
        
    end
end


%% Simulation Setup Parameter

SimSetup = struct();

SimSetup.Ts = 1e0; % Sampling time
SimSetup.K= 1e3; % Total number of time steps for simulation
SimSetup.T = SimSetup.Ts * SimSetup.K; % total simulation time
SimSetup.t  = (0:SimSetup.K-1)' * SimSetup.Ts; % Time vector

%% Simulate 'true' known input, u, and unknown input, d. 
u = zeros(1,SimSetup.K);

d = zeros(size(G,2),SimSetup.K); 
d(1,500:700) = 1*ones(size(d(1,500:700)));

for j = 100:800
    d(2,j) = 1/700*(j-100);
end
d(3,[500:549 600:649 700:749]) = 3*ones(size(d(3,[500:549 600:649 700:749])));
d(3,[550:599 650:699 750:799]) = -3*ones(size(d(3,[550:599 650:699 750:799])));
%% Generate the process noise and measurement noise vectors

SigmaQ = chol(Q);
w = (randn(SimSetup.K,size(A,1))*SigmaQ)';

SigmaR = chol(R);
v = (randn(SimSetup.K,size(C,1))*SigmaR)';

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
