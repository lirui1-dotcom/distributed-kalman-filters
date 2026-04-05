%% Initialization

clc;clear;close all;
% state space model--------------------------------------------------------
A = [0 -1 ; 1 0];
n = size (A,2);

H = eye(n);
m = size (A,1);

L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];

N=size(L,1);
beta=1;
L=beta*L;


% Covariance Matrices------------------------------------------------------

Q = diag([1,1]);
Ri=1;
R = diag([Ri,Ri,Ri,Ri,Ri,Ri,Ri,Ri]);     

% Initial Conditions ------------------------------------------------------
dt=0.01; % dt=0.001 is stable for a longer time, dt=0.02 is unstable when t→∞
t_f=20;
steps=round(t_f/dt);
x0 = [15 -10]'; %initial condition
Ad= (eye(n)+dt*A+dt^2/2*A^2+dt^3/6*A^3); % 3rd order discretization 
P_plus0 = [1 0 ;0 1];% initial covariance

% generate noise-----------------------------------------------------------
rng('default')
%rng('shuffle')
% W，n x k_steps
w=chol(Q)* randn(n,steps);
w=w*dt;
% V，n x k_steps
v = chol(R)* randn(n*N,steps);% shift v to right by 1


% Define matrices used to store all  data
% KF ---------------------------------------------------------------------- 
% x, used to store all states, ∈ R^(n x steps)
x = zeros(n,steps);

% % u, used to store inputs, ∈ R^(p x steps)
% u = zeros(p,steps);

% z, used to store measurements, ∈ R^(n x steps)
z = zeros(m,steps);
zh = zeros(m,steps);

% x_hat_plus，used to store posteriori estimates，∈ R^(n x steps)
x_hat_plus = zeros(n,steps);

% x_hat_minus_history，used to store priori estimates，∈ R^(n x steps)
x_hat_minus = zeros(n,steps);

% P_plus, used to store posteriori covariance matrices
P_plus = cell(1,steps); %cell array used to store all P^+

% P_plus, used to store priori covariance matrices
P_minus = cell(1,steps); %cell array used to store all P^-1



% True sates simulation 
x(:,1)=x0;
z(:,1)=H*x0+v(1); % z(1) is z[0]=Hx[0]+v[0]
for k = 2:steps
% state space euqation
   %x(:,k) = A * x(:,k-1) + B*u(:,k-1) +w(:,k-1);
   x(:,k) = Ad * x(:,k-1) + w(:,k-1); 
% measurement
   zh(:,k) = H * x(:,k); 
   z(:,k)= H * x(:,k)+v(1:2,k); 
end

%% Consensus on measurement
q = zeros(4*2,steps);
y = zeros(4*2,steps);
L = 1*L;
L_hat = kron(L,eye(m));

H1=[1 0];
H2=[0 1];

y(:,1)=[15
        0

        15
        0

        0
       -10
        
        0
       -10]';

u=[H1'*Ri^-1*zh(1,:); 
   H1'*Ri^-1*zh(1,:); 
   H2'*Ri^-1*zh(2,:);
   H2'*Ri^-1*zh(2,:)];


for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
   q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + (u(:,k-1)+v(:,k-1))  ) ; 
% y[k]=q[k]+z[k]
   y(:,k) = q(:,k) +  (u(:,k)+v(:,k)) ; 
end

% check disagreement-------------------------------------------------------
phi_y = zeros(1, steps); % Initialize phi to store results

for i = 1:steps
    yi = y(:, i);  % Extract the i-th column of y
    phi_y(i) = yi' * L_hat * yi;  % Compute phi for the i-th column
end

% plot---------------------------------------------------------------------
t_axis=(1:steps)*dt;
figure
hold on
plot(t_axis,y(1,:),'DisplayName','y1')  
%plot(y(2,:),'DisplayName','y2') 
plot(t_axis,y(3,:),'DisplayName','y3')
%plot(y(4,:),'DisplayName','y4')
plot(t_axis,y(5,:),'DisplayName','y5')
%plot(y(6,:),'DisplayName','y6')
plot(t_axis,y(7,:),'DisplayName','y7')
%plot(y(8,:),'DisplayName','y8')

plot(t_axis,zh(1,:)/2,'r--','DisplayName','z1 no measurement noise')
%plot(zh(2,:)+30,'r--','DisplayName','z2 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 

figure
plot(t_axis,phi_y)
title('Laplacian potential for y')

%% Consensus on covariance (only one col)
X = zeros(4*2,steps);
S = zeros(4*2,steps);
L = 1*L;
L_hat = kron(L,eye(m));

H1=[1 0];
H2=[0 1];

S(:,1)=[1
        0

        1
        0

        0
        1
        
        0
        1]';


Uc=[H1'*Ri^-1*H1; 
   H1'*Ri^-1*H1; 
   H2'*Ri^-1*H2;
   H2'*Ri^-1*H2];
U = Uc(:,1);

for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
  %X(:,k) = X(:,k-1)+ dt *(-L_hat)*(X(:,k-1) + U(:,k-1)  ) ; 
   X(:,k) = X(:,k-1)+ dt *(-L_hat)*(X(:,k-1) + U  ) ; 
% y[k]=q[k]+z[k]
   S(:,k) = X(:,k) + U ; 
end

%plot(z(1,:),'DisplayName','z1')
hold on
plot(S(1,:),'DisplayName','S1')
%plot(y(2,:),'DisplayName','S2')
plot(S(3,:),'DisplayName','S3')
%plot(y(4,:),'DisplayName','S4')
plot(S(5,:),'DisplayName','S5')
%plot(y(6,:),'DisplayName','S6')
plot(S(7,:),'DisplayName','S7')
%plot(y(8,:),'DisplayName','S8')

%plot(zh(1,:),'r--','DisplayName','z1 no measurement noise')
%plot(zh(2,:)+30,'r--','DisplayName','z2 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 
%% Consensus on covariance (both cols)
X = cell(1,steps);
S = cell(1,steps);
U = cell(1,steps); % we use cell in this case for convenience
L = 1*L;
L_hat = kron(L,eye(m));

H1=[1 0];
H2=[0 1];

Uc=[H1'*Ri^-1*H1; 
    H1'*Ri^-1*H1; 
    H2'*Ri^-1*H2;
    H2'*Ri^-1*H2];

for i=1:steps
    U{i}=Uc;
end

X{1}=zeros(8,1);
S{1}=Uc;

for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
  X{k} = X{k-1}+ dt *(-L_hat)*(X{k-1} + U{k-1}); 
   %X(:,k) = X(:,k-1)+ dt *(-L_hat)*(X(:,k-1) + U  ) ; 
% y[k]=q[k]+z[k]
   S{k} = X{k} + U{k} ; 
end

S_col_1 =zeros(8,steps);
S_col_2=zeros(8,steps);
% Loop through each cell in S
for i = 1:steps
    % Extract the columns of the matrix in each cell
    S_col_1(:,i) = S{i}(:, 1);
    S_col_2(:,i)=S{i}(:, 2);
end

% check disagreement-------------------------------------------------------
phi_S = zeros(2, steps); % Initialize phi to store results

for i = 1:steps  
    phi_S(1,i) = S_col_1(:,i)' * L_hat * S_col_1(:,i);  
    phi_S(2,i) = S_col_2(:,i)' * L_hat * S_col_2(:,i);  
end



t_axis=(1:steps)*dt;
figure
hold on
plot(t_axis,S_col_1(1,:),'DisplayName','S1(:,1)')
plot(t_axis,S_col_2(1,:),'DisplayName','S1(:,2)')

plot(t_axis,S_col_1(3,:),'DisplayName','S3(:,1)')
plot(t_axis,S_col_2(3,:),'DisplayName','S3(:,2)')

plot(t_axis,S_col_1(5,:),'DisplayName','S5(:,1)')
plot(t_axis,S_col_2(5,:),'DisplayName','S5(:,2)')

plot(t_axis,S_col_1(7,:),'DisplayName','S7(:,1)')
plot(t_axis,S_col_2(7,:),'DisplayName','S7(:,2)')

%plot(zh(1,:),'r--','DisplayName','z1 no measurement noise')
%plot(zh(2,:)+30,'r--','DisplayName','z2 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 

figure
plot(t_axis,phi_S(1,:),'DisplayName','S1(:,1) disagreement')
hold on
plot(t_axis,phi_S(2,:),'DisplayName','S2(:,1) disagreement')
legend 

