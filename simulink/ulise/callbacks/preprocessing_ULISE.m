function [B, D, input_data] = preprocessing_ULISE(A, B, C, D, G, K, Ts, u, d, filename)
%PREPROCESSING_ULISE Prepare LTI system inputs for ULISE Simulink model
%
% This function creates timeseries objects for known (u) and unknown (d)
% inputs, handles zero-input edge cases, and prepares a Simulink dataset
% compatible with Signal Editor blocks. The dataset is saved to
% a MAT-file for reuse.
%
% Inputs:
%   A       - State transition matrix (n x n)
%   B       - Control input matrix (n x m)
%   C       - Measurement matrix (l x n)
%   D       - Direct feedthrough from u (l x m)
%   G       - Disturbance input matrix (n x p)
%   K       - Total number of simulation time steps (scalar)
%   Ts      - Sampling time (scalar)
%   u       - Known input matrix (m x K)
%   d       - Unknown input matrix (p x K)
%   filename - String specifying output MAT-file name
%
% Outputs:
%   B          - Possibly updated control input matrix (n x m or n x 1 if m=0)
%   D          - Possibly updated feedthrough matrix (l x m or l x 1 if m=0)
%   input_data - Simulink.SimulationData.Dataset containing u and d timeseries
%
%   Authorship:
%       D. Lawson, S. Macero, E. Gah, and S.Z. Yong (2025). SISE Toolbox
%       for MATLAB (https://github.com/dalawson0/SISE.git), GitHub. 

arguments (Input)
    A
    B
    C
    D
    G
    K
    Ts
    u
    d
    filename
end


%% Make timeseries input data


n = size(A,1); % Number of states
m = size(B,2); % Number of known inputs
l = size(C,1); % Number of Outputs
p = size(G,2); % Number of unknown inputs

t  = (0:K-1)' * Ts; % Time vector

% Replaces the input (and B and D matrices) with scalar dummy zero input
if m == 0
    u = zeros(1, K);
    B = zeros(n, 1);
    D = zeros(l, 1);
end

u_ts = timeseries(u, t);
d_ts = timeseries(d, t);

%% Create dataset to utilize Signal Editor block

input_data = Simulink.SimulationData.Dataset;
input_data = input_data.addElement(setinterpmethod(d_ts,'zoh'), 'd');
input_data = input_data.addElement(setinterpmethod(u_ts,'zoh'), 'u');


% Save dataset to file
full_path = fullfile('./Examples/Simulink/', filename);
addpath('./Examples/Simulink/');
fprintf('Saving to: %s \n', full_path);
save(full_path, 'input_data');
end