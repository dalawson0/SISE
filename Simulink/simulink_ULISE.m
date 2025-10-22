% Daniel Lawson
% Oct, 10, 2025
% Script for Simulink implementation of ULISE

% Define state-space matrices
A = [0.5    2   0   0   0;
     0      0.2 1   0   1; 
     0      0   0.3 0   1; 
     0      0   0   0.7 1;
     0      0   0   0   0.1]; 

B = zeros(size(A,1),1);

C = blkdiag(1,1,1,1,1);

D = zeros(size(C,1),1);

G = [1  0   -0.3;
     1   0   0;
     0   0   0;
     0   0   0;
     0   0   0];

H = [0  0   1;
     0  0   0; 
     0  1   0; 
     0  0   0; 
     0  0   0]; 

R = 1e-2* [1    0   0   0.5 0;
           0    1   0   0   0.3; 
           0    0   1   0   0; 
           0.5  0   0   1   0; 
           0    0.3 0   0   1];

Q = 1e-4* [1    0   0   0   0; 
           0    1   0.5 0   0; 
           0    0.5 1   0   0; 
           0    0   0   1   0; 
           0    0   0   0   1]; 

% K = 1000;

x0 = zeros(size(A,1),1);        % Initial true state 
% xhat0 = zeros(size(A,1),1);     % Initial state estimate 
% P_x0 = eye(size(A,1))*1e6;      % Initial state error covariance 

K = 1000; % Total time steps

% Next, we define the known input  and the true disturbance trajectory  as column vectors for all time steps. 
u = zeros(1,K);

d = zeros(size(G,2),K); 
d(1,500:700) = 1*ones(size(d(1,500:700)));

for j = 100:800
    d(2,j) = 1/700*(j-100);
end
d(3,[500:549 600:649 700:749]) = 3*ones(size(d(3,[500:549 600:649 700:749])));
d(3,[550:599 650:699 750:799]) = -3*ones(size(d(3,[550:599 650:699 750:799])));

Ts = 1e0; % sample time
t = (0:K - 1)' * Ts; % create time vector 
d_timeseries = timeseries(d, t'); % create 4x1000
u_timeseries = timeseries(u, t'); % create 2x1000

% The process and measurement noise covariance matrices  and  are used to generate the process noise  and measurement noise . 
SigmaQ = chol(Q);
w = (randn(K,size(A,1))*SigmaQ)';

SigmaR = chol(R);
v = (randn(K,size(C,1))*SigmaR)';

% 
w_timeseries = timeseries(w, t'); % create 6x1000
v_timeseries = timeseries(v, t'); % create 6x1000

% time simulation
total_time = Ts * K;

x = zeros(size(A,1),K);
y = zeros(size(C,1),K);

x(:,1)=x0;

for i=1:K
    y(:,i)=C*x(:,i)+D*u(i)+H*d(:,i)+v(:,i);
    if i<K
        x(:,i+1)=A*x(:,i)+B*u(i)+G*d(:,i)+w(:,i);
    end
end

x(:,end)