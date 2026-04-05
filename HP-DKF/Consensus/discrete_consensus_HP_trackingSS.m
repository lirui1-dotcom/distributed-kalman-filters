%% case 1:each node track full states
clc;clear;close all;
% state space model--------------------------------------------------------
A = [0 -1;1 0];
n = size (A,2);

H = eye(n);
m = size (A,1);

L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];
N=size(L,1);
% Covariance Matrices------------------------------------------------------

Q = diag([1,1]);
Ri=1;
R = diag([Ri,Ri,Ri,Ri,Ri,Ri,Ri,Ri]);     

% Initial Conditions ------------------------------------------------------
dt=0.02; % dt=0.001 is stable for a longer time, dt=0.02 is unstable when t→∞
t_f=20;
steps=round(t_f/dt);
x0 = [15 -10]'; %initial condition
Ad= (dt*A+eye(n)); % discretization using euler's method

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

% K, used to store all kalman gains
K = cell(1,steps);

% χ2 ----------------------------------------------------------------------
% V, used to store all inovation vectors
Nu = zeros(m,steps);
% S, used to store all inovation covariances matrices
S = cell(1,steps);
% Chi-square , used to store all χ2 
Chi = zeros(1,steps);

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


q = zeros(4*2,steps);
y = zeros(4*2,steps);
L = 1*L;
L_hat = kron(L,eye(m));


y(:,1)=[15+10 -10+10 15+20 -10+20 15+30 -10+30  15+60  -10+60]';
u=[zh(1,:)+10;  zh(2,:)+10; zh(1,:)+20 ; zh(2,:)+20 ; zh(1,:)+30 ; zh(2,:)+30 ; zh(1,:)+60 ; zh(2,:)+60];

for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
   q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + (u(:,k-1)+v(:,k-1))  ) ; 
% y[k]=q[k]+z[k]
   y(:,k) = q(:,k) +  (u(:,k)+v(:,k)) ; %+ [v1(:,k); v2(:,k); v3(:,k); v4(:,k)];
end

%plot(z(1,:),'DisplayName','z1')
hold on
plot(y(1,:),'DisplayName','y1')
plot(y(2,:),'DisplayName','y2')
plot(y(3,:),'DisplayName','y3')
plot(y(4,:),'DisplayName','y4')
plot(y(5,:),'DisplayName','y5')
plot(y(6,:),'DisplayName','y6')
plot(y(7,:),'DisplayName','y7')
plot(y(8,:),'DisplayName','y8')

plot(zh(1,:)+30,'r--','DisplayName','z1 no measurement noise')
plot(zh(2,:)+30,'r--','DisplayName','z2 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 

%% case 2:each node track partial states
clc;clear;close all;
% state space model--------------------------------------------------------
A = [0 -1;1 0];
n = size (A,2);

H = eye(n);
m = size (A,1);

L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];
beta=10;
L=beta*L;
N=size(L,1);
% Covariance Matrices------------------------------------------------------

Q = diag([1,1]);
Ri=1;
R = diag([Ri,Ri,Ri,Ri,Ri,Ri,Ri,Ri]);     

% Initial Conditions ------------------------------------------------------
dt=0.001; % dt=0.001 is stable for a longer time, dt=0.02 is unstable when t→∞
t_f=20;
steps=round(t_f/dt);
x0 = [15 -10]'; %initial condition
Ad= (eye(n)+dt*A+dt^2/2*A^2+dt^3/6*A^3); % 3rd order discretization 

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

% True sates simulation 

x(:,1)=x0;
z(:,1)=H*x0+v(1); % z(1) is z[0]=Hx[0]+v[0]
for k = 2:steps
% state space euqation
   %x(:,k) = A * x(:,k-1) + B*u(:,k-1) +w(:,k-1);
   x(:,k) = Ad * x(:,k-1) + w(:,k-1); 
% measurement
   zh(:,k) = H * x(:,k); 
   z(:,k)  = H * x(:,k)+v(1:2,k); 
end


q = zeros(4*2,steps);
y = zeros(4*2,steps);
L = 1*L;
L_hat = kron(L,eye(m));


y(:,1)=[15+10 
       %-10+10 
        0

        15+30 
       %-10+20
        0

        0
       % 15+30 
       -10+30
        
        0
        %15+60  
       -10+60]';
u=[zh(1,:)+10 ; 
%  zh(2,:)+10 ; 
   zeros(1,steps);

   zh(1,:)+30 ; 
%  zh(2,:)+20 ;
   zeros(1,steps);

   zeros(1,steps);
%  zh(1,:)+30 ;
   zh(2,:)+30 ;
   
   zeros(1,steps);
%  zh(1,:)+60 ;
   zh(2,:)+60 ];

for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
   %q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + (u(:,k-1)+v(:,k-1))  ) ; 
   q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + (u(:,k-1)) ); 
% y[k]=q[k]+z[k]
   y(:,k) = q(:,k) +  (u(:,k)) ; 
   %y(:,k) = q(:,k) +  (u(:,k)+v(:,k)) ; 
end

%plot(z(1,:),'DisplayName','z1')
hold on
plot(y(1,:),'DisplayName','y1')
%plot(y(2,:),'DisplayName','y2')
plot(y(3,:),'DisplayName','y3')
%plot(y(4,:),'DisplayName','y4')
plot(y(5,:),'DisplayName','y5')
%plot(y(6,:),'DisplayName','y6')
plot(y(7,:),'DisplayName','y7')
%plot(y(8,:),'DisplayName','y8')

plot(0.5*zh(1,:)+10,'r--','DisplayName','z1 no measurement noise')
%plot(zh(2,:)+30,'r--','DisplayName','z2 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 

%% case 3:each node track Hi'*Ri^-1*Zi
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

% note that this case each node can observe all states, then each node reach H'*R^-1*Z, where H1=H2=...=Hn=H
% in this case because each node has similar values, beta=0.01 also works
% but this is not ideal for distributed, in distributed case some nodes can only measure some parts of the states
% i.g. u=[z1 0 z1 0 0 z2 0 z2]', node 1 3 5 7 will try to reach consensus
% with y1=y3=y5=y7= (z1+z1)/4=(z1+z2)/2, but node 5 and node 7 are not meausring z1, y3[0]=y7[0]=0, which require fast convergence rate

% u=[H'*(diag([Ri,Ri]))^-1*zh; 
%    H'*(diag([Ri,Ri]))^-1*zh; 
%    H'*(diag([Ri,Ri]))^-1*zh;
%    H'*(diag([Ri,Ri]))^-1*zh];

for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
   q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + (u(:,k-1)+v(:,k-1))  ) ; 
   %q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + (u(:,k-1)) ); 
% y[k]=q[k]+z[k]
   %y(:,k) = q(:,k) +  (u(:,k)) ; 
   y(:,k) = q(:,k) +  (u(:,k)+v(:,k)) ; 
end

%plot(z(1,:),'DisplayName','z1')
hold on
plot(y(1,:),'DisplayName','y1') 
%plot(y(2,:),'DisplayName','y2')
plot(y(3,:),'DisplayName','y3')
%plot(y(4,:),'DisplayName','y4')
plot(y(5,:),'DisplayName','y5')
%plot(y(6,:),'DisplayName','y6')
plot(y(7,:),'DisplayName','y7')
%plot(y(8,:),'DisplayName','y8')

plot(zh(1,:)/2,'r--','DisplayName','z1 no measurement noise')
%plot(zh(2,:)+30,'r--','DisplayName','z2 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 
