function [xhat_UL,Px_UL,dhat_UL,Pd_UL,Pxd_UL] = ULISE(A,B,C,D,G,H,Q,R,K,u,y,xhat0,P_x0)
% ULISE Unified Linear Input & State Estimator for Linear Time-Invariant
%       Systems (Zero-Step State, One-Step Input Delay) 
%
%   Syntax:
%       [xhat_UL, Px_UL, dhat_UL, Pd_UL, Pxd_UL] = ULISE(A,B,C,D,G,H,Q,R,K,u,y,xhat0,P_x0)
%
%   Description:
%       ULISE computes minimum-variance, unbiased estimates of the state
%       x_k and unknown input (disturbance) d_{k-1} with one-step delay,
%       for a discrete-time linear time-invariant system: 
%
%           x_{k+1} = A x_k + B u_k + G d_k + w_k   {State equation}
%           y_k     = C x_k + D u_k + H d_k + v_k   {Measurement equation}
%
%       where u_k is the known input vector, w_k is the zero-mean process
%       noise vector with covariance Q, and v_k is zero-mean measurement
%       noise vector with covariance R.
%
%       The state and unknown input estimates are computed using the
%       equations presented in Section 4 of [1]. These include the unknown
%       input estimation, state prediction / time update, and state measurement update
%       steps. For complete derivations and definitions of intermediate
%       matrices (e.g., M_{1,k}, M_{2,k}, and L_k — implemented in code as
%       M1, M_UL, and L_UL), refer to [1].
%
%   Inputs:
%       A           : State transition matrix (n-by-n)
%       B           : Control input matrix (n-by-m)
%       C           : Measurement matrix (l-by-n)
%       D           : Measurement feedthrough from u_k (l-by-m)
%       G           : Disturbance input matrix (n-by-p)
%       H           : Disturbance feedthrough to output (l-by-p)
%       Q           : Process noise covariance (n-by-n)
%       R           : Measurement noise covariance (l-by-l)
%       K           : Number of time steps (scalar)
%       u           : Known control inputs; each column is u_k at time k (m-by-K)
%       y           : Measured outputs; each column is y_k at time k (l-by-K)
%       xhat0       : Initial state estimate (n-by-1)
%       P_x0        : Initial state error covariance (n-by-n)
%
%   Outputs:
%       xhat_UL     : Estimated states; each column is xhat_{k|k} at each time step k (n-by-K)
%       Px_UL       : State error covariance; block entries correspond to P^x_{k|k} at each time step k (n-by-n-by-K)
%       dhat_UL     : Estimated one-step delay unknown input; each column is dhat_{k-1} at each time step k-1 (p-by-p-by-(K-1))
%       Pd_UL       : Unknown input error covariance; block entries correspond to P^d_{k-1} at each time step k-1 (p-by-p-by-(K-1))
%       Pxd_UL      : Cross-covariance between state and input estimates;block entries correspond to P^{xd}_{k-1} at each time step k-1 (n-by-p-by-(K-1))
%
%   Notes:
%       - Unobservable components of d_k (denoted d_{2,k}) are estimated one
%         time step later if rank(H) ≠ full rank.
%
%   References:
%       [1] Yong, S. Z., Zhu, M., & Frazzoli, E. (2016). A unified filter
%       for simultaneous input and state estimation of linear discrete-time
%       stochastic systems. Automatica, 63, 321–329.
%       https://doi.org/10.1016/j.automatica.2015.10.040 
%       (Also available as a preprint at https://arxiv.org/abs/1309.6627)
%
%   Authorship:
%       D. Lawson, S. Macero, E. Gah, and S.Z. Yong (2025). SISE Toolbox
%       for MATLAB (https://github.com/dalawson0/SISE.git), GitHub. 
%   
%   Version History:
%       Introduced in SISE v1.0.0

n=size(A,1);
l=size(C,1);
p=size(G,2);
r=rank(H);

% Generate time steps
k=1:K;

% Coordinate transformation
[U,S,V]=svd(H);
U1=U(:,1:r);
U2=U(:,r+1:end);
Sigma=S(1:r,1:r);
V1=V(:,1:r);
V2=V(:,r+1:end);
T1=[eye(r) -U1'*R*U2*pinv(U2'*R*U2)]*[U1';U2'];
T2=U2';
C1=T1*C;
C2=T2*C;
D1=T1*D;
D2=T2*D;
G1=G*V1;
G2=G*V2;
M1=pinv(Sigma);
R1=T1*R*T1';
R2=T2*R*T2';
Ahat=A-G1*M1*C1;
Qhat=G1*M1*R1*M1'*G1'+Q;

% Check Assumptions (i) and (ii)
if license("test", "Control_Toolbox")
    t_zeros = tzero(A,G,C,H)
    if rank([A G;C H]) ~= n+p
        error('Error. System is not strongly detectable.')
    end
    for i = 1:length(t_zeros)
        if abs(t_zeros(i))>1
            error('Error. System is not strongly detectable.')
        end
    end
else
    warning("Control System Toolbox not available. Condition (i) is not checked.")
end

if rank(C2*G2)<p-r
    error('Error. Delay greater or equal to 1. See: Yong, S.Z., Zhu, M., and Frazzoli, E. (2015). Simultaneous input and state estimation with a delay. IEEE Conference on Decision and Control, Osaka, Japan, pp. 468-475')
end

% Output decoupling
z1=zeros(r,K);
z2=zeros(l-r,K);
for i = 1:K
    z1(:,i)=T1*y(:,i);
    z2(:,i)=T2*y(:,i);
end

%% Initialization
xhat_UL=zeros(n,length(k));
xhat_p_UL=zeros(n,length(k));
xhat_star_UL=zeros(n,length(k));

dhat_UL=zeros(p,length(k));
Px_UL=zeros(n,n,length(k));
Pxd_UL=zeros(n,p,length(k));
Pd_UL=zeros(p,p,length(k));
Px_star_UL=zeros(n,n,length(k));
P_tilde_UL=zeros(n,n,length(k));

R1_tilde_UL=zeros(r,r,length(k));
R2_tilde_UL=zeros(l-r,l-r,length(k));
R2_tilde_star_UL=zeros(l-r,l-r,length(k));

Pd1_UL=zeros(r,r,length(k));
Pd2_UL=zeros(p-r,p-r,length(k));
Pd12_UL=zeros(r,p-r,length(k));
Pxd1_UL=zeros(n,r,length(k));
Pxd2_UL=zeros(n,p-r,length(k));

M_UL=zeros(p-r,l-r,length(k));
L_UL=zeros(n,l-r,length(k));

dhat1_UL=zeros(r,length(k));
dhat2_UL=zeros(p-r,length(k));

xhat_UL(:,1) = xhat0;
Px_UL(:,:,1) = P_x0;
dhat1_UL(:,1)=M1*(z1(:,1)-C1*xhat_UL(:,1)-D1*u(:,1));
Pd1_UL(:,:,1) = M1*(C1*Px_UL(:,:,1)*C1'+R1)*M1;
Pxd1_UL(:,:,1) = -Px_UL(:,:,1)*C1'*M1';

% Filter dynamics (ULISE)
for i=k(1:end-1)
    % Estimation of d2 and d
    P_tilde_UL(:,:,i+1)=Ahat*Px_UL(:,:,i)*Ahat'+Qhat; 
    R2_tilde_UL(:,:,i+1)=C2*P_tilde_UL(:,:,i+1)*C2'+R2;
    Pd2_UL(:,:,i)=inv(G2'*C2'*inv(R2_tilde_UL(:,:,i+1))*C2*G2);
    M_UL(:,:,i+1)=Pd2_UL(:,:,i)*G2'*C2'*inv(R2_tilde_UL(:,:,i+1));
    xhat_p_UL(:,i+1)=A*xhat_UL(:,i)+B*u(:,i)+G1*dhat1_UL(:,i);
    dhat2_UL(:,i)=M_UL(:,:,i+1)*(z2(:,i+1)-C2*xhat_p_UL(:,i+1)-D2*u(:,i+1));
    dhat_UL(:,i)=V*[dhat1_UL(:,i);dhat2_UL(:,i)];
    Pd12_UL(:,:,i)=-Pxd1_UL(:,:,i)'*A'*C2'*M_UL(:,:,i+1)'-Pd1_UL(:,:,i)*G1'*C2'*M_UL(:,:,i+1)';
    Pxd2_UL(:,:,i)=-Px_UL(:,:,i)*A'*C2'*M_UL(:,:,i+1)'-Pxd1_UL(:,:,i)*G1'*C2'*M_UL(:,:,i+1)';
    
    % State prediction / Time update 
    xhat_star_UL(:,i+1)=xhat_p_UL(:,i+1)+G2*dhat2_UL(:,i);
    Px_star_UL(:,:,i+1)=G2*M_UL(:,:,i+1)*R2*M_UL(:,:,i+1)'*G2'+(eye(n)-G2*M_UL(:,:,i+1)*C2)*P_tilde_UL(:,:,i+1)*(eye(n)-G2*M_UL(:,:,i+1)*C2)'; 
    R2_tilde_star_UL(:,:,i+1)=C2*Px_star_UL(:,:,i+1)*C2'+R2-C2*G2*M_UL(:,:,i+1)*R2-R2*M_UL(:,:,i+1)'*G2'*C2';
    
    % Measurement update
    R_pinv=pinv(R2_tilde_star_UL(:,:,i+1));
    L_UL(:,:,i+1)=(Px_star_UL(:,:,i+1)*C2'-G2*M_UL(:,:,i+1)*R2)*R_pinv;
    xhat_UL(:,i+1)=xhat_star_UL(:,i+1)+L_UL(:,:,i+1)*(z2(:,i+1)-D2*u(i+1)-C2*xhat_star_UL(:,i+1));
    Px_UL(:,:,i+1)=(eye(n)-L_UL(:,:,i+1)*C2)*Px_star_UL(:,:,i+1)*(eye(n)-L_UL(:,:,i+1)*C2)'+L_UL(:,:,i+1)*R2*L_UL(:,:,i+1)'+L_UL(:,:,i+1)*R2*M_UL(:,:,i+1)'*G2'*(eye(n)-L_UL(:,:,i+1)*C2)'+(eye(n)-L_UL(:,:,i+1)*C2)*G2*M_UL(:,:,i+1)*R2*L_UL(:,:,i+1)';
    
    % Estimation of d1
    R1_tilde_UL(:,:,i+1)=C1*Px_UL(:,:,i+1)*C1'+R1;
    Pd1_UL(:,:,i+1)=M1*R1_tilde_UL(:,:,i+1)*M1;
    dhat1_UL(:,i+1)=M1*(z1(:,i+1)-C1*xhat_UL(:,i+1)-D1*u(:,i+1));
    Pxd1_UL(:,:,i+1)=-Px_UL(:,:,i+1)*C1'*M1';
    Pd_UL(:,:,i)=V*[Pd1_UL(:,:,i) Pd12_UL(:,:,i); Pd12_UL(:,:,i)' Pd2_UL(:,:,i)]*V';
    
    % Estimation of Pxd
    Pxd_UL(:,:,i)=Pxd1_UL(:,:,i)*V1'+Pxd2_UL(:,:,i)*V2';
end 
