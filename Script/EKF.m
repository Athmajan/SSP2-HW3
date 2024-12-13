%% Lab 3 - Extended Kalman Filter
close all;
clear;
clc;

rng(1);

N = 100;

delta = 1;
%%-----------------------------------
% % Case I
% sig_u_sqr = 0.0001;
% sig_r_sqr = 0.1;
% sig_b_sqr = 0.01;

% %%-----------------------------------
% % % Case II
% sig_u_sqr = 0.01;
% sig_r_sqr = 0.1;
% sig_b_sqr = 0.01;
% 
% 
% %%-----------------------------------
% % %Case III
sig_u_sqr = 0.0001;
sig_r_sqr = 1;
sig_b_sqr = 0.5;

A = eye(4) + diag([delta delta],2);
auxvar = [0 0 sig_u_sqr sig_u_sqr];
Q = diag(auxvar);
C = diag([sig_r_sqr sig_b_sqr]);

Un = sqrt(sig_u_sqr)*randn(2,N);
U = [zeros(2,N); Un];
R_ideal = zeros(2,N);

vn =[-0.2; 0.2]; % vn[-1]
rn = [10; -5];  % rn[-1]

rn_ideal = rn;
R_ideal(:,1) = rn;

% % Change Initial State
%%----------------------
S(:,1)= [rn; vn]; % Given in the task
%S(:,1)= [0;0;-0.1;0.1]; % My change to init state

h_sn(:,1) = [ sqrt(S(1,1)^2 + S(2,1)^2); atan( S(2,1)/S(1,1) ) ];
for n=2:N
    rn_ideal = rn_ideal + [-0.2; 0.2]*delta; 
    R_ideal(:,n) = rn_ideal;
    S(:,n) = A*S(:,n-1) + U(:,n);
    h_sn(:,n) = [ sqrt(S(1,n)^2 + S(2,n)^2); atan2( S(2,n),S(1,n) ) ];
end

w = [sqrt(sig_r_sqr)*randn(1,N); sqrt(sig_b_sqr)*randn(1,N)];
xn = h_sn + w;

R_obs = [xn(1,:).*cos(xn(2,:)); xn(1,:).*sin(xn(2,:))];

figure
plot(R_ideal(1,:),R_ideal(2,:),'--k',S(1,:),S(2,:),'-k',R_obs(1,:),R_obs(2,:),'-b');
xlabel('$r_x$','interpreter','latex')
ylabel('$r_y$','interpreter','latex')
legend('Ideal track', 'True track', 'Observed track')
grid on;

figure
plot(0:N-1,h_sn(1,:)); 
xlabel('n','interpreter','latex')
ylabel('$R[n]$','interpreter','latex')
title('Range','interpreter','latex')
grid on;

figure
plot(0:N-1,h_sn(2,:)*180/pi); 
xlabel('n','interpreter','latex')
ylabel('$\beta [n]$ (degrees)','interpreter','latex')
title('Bearing','interpreter','latex')
grid on;
%% Extended Kalman filter

S_hat = zeros(4,N);

% Initial Estimates to change
%S_hat(:,1) = [5 5 0 0]';
S_hat(:,1) = [10 10 0.1 0.1]';

M = zeros(4,4,N);
M(:,:,1) = eye(4);

for n=2:N
    S_hat(:,n) = A*S_hat(:,n-1);
    M(:,:,n) = A*M(:,:,n-1)*A' + Q;
    
    rxn = S_hat(1,n);
    ryn = S_hat(2,n);
    lenRxy = sqrt(rxn^2 + ryn^2);

    H = [rxn/lenRxy ryn/lenRxy 0 0;
         -1*ryn/lenRxy^2 rxn/lenRxy^2 0 0];

    K_num = M(:,:,n)*H.';
    K_den = C + H*M(:,:,n)*H.';
    K = K_num/K_den;
    
    hsn = [lenRxy; atan2(ryn,rxn)];
    S_hat(:,n) = S_hat(:,n) + K*( xn(:,n) - hsn );


    M(:,:,n) = (eye(4)-K*H)*M(:,:,n);

end


figure
plot(S(1,:),S(2,:),'k--',S_hat(1,:),S_hat(2,:),'k','LineWidth',2);
xlabel(['Kalman estimate and true ' 'r_{x}[n]'])
ylabel(['Kalman estimate and true ' 'r_{y}[n]'])
legend('True track','EKF estimate')

grid on
figure
plot(0:N-1,squeeze(M(1,1,:)),'k');
xlabel('Sample Number n')
ylabel('Minimum MSE for r_{x}[n]')



grid on
figure
plot(0:N-1,squeeze(M(2,2,:)),'k');
xlabel('Sample Number n')
ylabel('Minimum MSE for r_{y}[n]')














