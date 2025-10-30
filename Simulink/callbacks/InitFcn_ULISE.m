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

% Ahat = A - G1*M1*C1;
% Qhat = G1*M1*R1*M1'*G1' + Q;

T1 = U1' - U1'*R*U2*inv(U2'*R*U2)*U2';
T2 = U2';

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

% if rank(C2*G2)??
PlantParams = struct();

PlantParams.A = Simulink.Parameter;
PlantParams.A.Value = [ 0.5    2   0   0   0;
                        0      0.2 1   0   1; 
                        0      0   0.3 0   1; 
                        0      0   0   0.7 1;
                        0      0   0   0   0.1 ]; 
PlantParams.A.Description = "State Transition Matrix";

PlantParams.B = Simulink.Parameter;
PlantParams.B.Value = zeros(size(PlantParams.A.Value,1),1);

PlantParams.C = Simulink.Parameter;
PlantParams.C.Value = blkdiag(1,1,1,1,1);

PlantParams.D = Simulink.Parameter; 
PlantParams.D.Value = zeros(size(PlantParams.C.Value,1),1);

PlantParams.G = Simulink.Parameter;
PlantParams.G.Value = [ 1  0   -0.3;
                        1   0   0;
                        0   0   0;
                        0   0   0;
                        0   0   0] ;
p = size(PlantParams.G.Value, 2);

PlantParams.H = Simulink.Parameter;
PlantParams.H.Value = [ 0  0   1;
                        0  0   0; 
                        0  1   0; 
                        0  0   0; 
                        0  0   0 ]; 

PlantParams.x0 = Simulink.Parameter;
PlantParams.x0.Value = zeros(size(PlantParams.A.Value,1),1);
PlantParams.x0.Description = "Initial true state";

% % Initialize the PlantParams structure with the defined Parameter
% assignin('base', 'PlantParams', PlantParams);
%% Define ULISE Parameter

ULISEParams = struct();

ULISEParams.R = Simulink.Parameter;
ULISEParams.R.Value = 1e-2 * [  1    0   0   0.5 0;
                                0    1   0   0   0.3; 
                                0    0   1   0   0; 
                                0.5  0   0   1   0; 
                                0    0.3 0   0   1  ];
ULISEParams.R.Description = "Measurement Noise Error Covariance";


ULISEParams.Q = Simulink.Parameter;
ULISEParams.Q.Value = 1e-4 * [  1    0   0   0   0; 
                                0    1   0.5 0   0; 
                                0    0.5 1   0   0; 
                                0    0   0   1   0; 
                                0    0   0   0   1  ]; 
ULISEParams.Q.Description = "Process Noise Error Covariance";

% Singular Value Decomposition of H
r = rank(PlantParams.H.Value);
[U, S, V] = svd(PlantParams.H.Value);

ULISEParams.U1 = Simulink.Parameter; 
ULISEParams.U1.Value = U(:, 1:r);

ULISEParams.U2 = Simulink.Parameter; 
ULISEParams.U2.Value = U(:, r+1:end);

ULISEParams.Sigma = Simulink.Parameter; 
ULISEParams.Sigma.Value = S(1:r, 1:r);

ULISEParams.V1 = Simulink.Parameter; 
ULISEParams.V1.Value = V(:, 1:r);

ULISEParams.V2 = Simulink.Parameter; 
ULISEParams.V2.Value = V(:, r+1:end);

ULISEParams.xhat0 = Simulink.Parameter;
ULISEParams.xhat0.Value = zeros(size(PlantParams.A.Value,1),1);
ULISEParams.xhat0.Description = "Initial state estimate";

ULISEParams.P_x0 = Simulink.Parameter;
ULISEParams.P_x0.Value = eye(size(PlantParams.A.Value,1))*1e6;
ULISEParams.P_x0.Description = "Initial state error covariance";

%% Simulation Setup Parameter

SimSetup = struct();

SimSetup.Ts = 1e0; % Time step
SimSetup.K= 1e3; % Total number of time steps for simulation
SimSetup.T = SimSetup.Ts * SimSetup.K; % total simulation time
SimSetup.t  = (0:SimSetup.K-1)' * SimSetup.Ts; % Time vector
% assignin('base', 'SimSetup', SimSetup); % Assign SimSetup to the base workspace

%% Simulate 'true' known input, u, and unknown input, d. 
u = zeros(1,SimSetup.K);

d = zeros(size(PlantParams.G.Value,2),SimSetup.K); 
d(1,500:700) = 1*ones(size(d(1,500:700)));

for j = 100:800
    d(2,j) = 1/700*(j-100);
end
d(3,[500:549 600:649 700:749]) = 3*ones(size(d(3,[500:549 600:649 700:749])));
d(3,[550:599 650:699 750:799]) = -3*ones(size(d(3,[550:599 650:699 750:799])));
%% Generate the process noise and measurement noise vectors

SigmaQ = chol(ULISEParams.Q.Value);
w = (randn(SimSetup.K,size(PlantParams.A.Value,1))*SigmaQ)';

SigmaR = chol(ULISEParams.R.Value);
v = (randn(SimSetup.K,size(PlantParams.C.Value,1))*SigmaR)';

%% Make timeseries data 

InputData = struct();

InputData.u=setinterpmethod(timeseries(u, SimSetup.t),'zoh');
% InputData.u = timeseries(u, SimSetup.t', "");
InputData.u.Name = "timeseries known input";
% u_stairs = setinterpmethod(, 'zoh');

% InputData.d = timeseries(d, SimSetup.t');
InputData.d=setinterpmethod(timeseries(d, SimSetup.t),'zoh');
InputData.d.Name = "timeseries unknown input";

% Add Noise timeseries to Plant Struct
PlantParams.w = timeseries(w, SimSetup.t');
PlantParams.w.Name = "timeseries process noise";

PlantParams.v = timeseries(v, SimSetup.t');
PlantParams.v.Name = "timeseries measurement noise";
