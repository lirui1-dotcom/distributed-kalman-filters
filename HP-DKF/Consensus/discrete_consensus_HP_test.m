%% case 1: z=[sin+10 sin+20 sin+30 sin+40]' 
clc;clear;close all
dt=0.02; 
t_f=20;
steps=round(t_f/dt);

L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];
q = zeros(4,steps);
y = zeros(4,steps);
L = 1*L;
%L_hat = kron(L,eye(m));


y(:,1)=[10 20 30 60]';
z=sin(10*dt*(1:steps));

%y(:,1) = q(:,1)+[z(:,1) z(:,1) z(:,1) z(:,1)]'; 

for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
  % q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + [  z(:,k-1)+10 ; z(:,k-1)+20 ;z(:,k-1)+30 ;z(:,k-1)+60]) ;
    q(:,k) = q(:,k-1)+ dt *(-L)*(q(:,k-1) + [  z(k-1)+10 ; z(k-1)+20 ;z(k-1)+30 ;z(k-1)+60]) ;
% y[k]=q[k]+z[k]
   y(:,k) = q(:,k) + [z(k)+10; z(k)+20; z(k)+30; z(k)+60]; %+ [v1(:,k); v2(:,k); v3(:,k); v4(:,k)];
end

%plot(z(1,:),'DisplayName','z1')
hold on
plot(y(1,:),'DisplayName','y1')
plot(y(2,:),'DisplayName','y2')
plot(y(3,:),'DisplayName','y3')
plot(y(4,:),'DisplayName','y4')
plot(z+30,'DisplayName','z1 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 
ylim([0,50])
xlim([0 200])
%% case 2: z=[sin+10+v1 sin+20+v2 sin+30+v3 sin+40+v4]'
clc;clear;close all
dt=0.02; 
t_f=20;
steps=round(t_f/dt);

L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];
q = zeros(4,steps);
y = zeros(4,steps);
L = 1*L;
%L_hat = kron(L,eye(m));


y(:,1)=[10 20 30 60]';
z=sin(10*dt*(1:steps));

%generate noise
R = diag([0.1,0.1,0.1,0.1]);      
v = chol(R)* randn(4,steps);

%y(:,1) = q(:,1)+[z(:,1) z(:,1) z(:,1) z(:,1)]'; 

for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
  % q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + [  z(:,k-1)+10 ; z(:,k-1)+20 ;z(:,k-1)+30 ;z(:,k-1)+60]) ;
    q(:,k) = q(:,k-1)+ dt *(-L)*(q(:,k-1) + [z(k-1)+10+v(1,k) ; z(k-1)+20+v(2,k) ;z(k-1)+30+v(3,k) ;z(k-1)+60+v(4,k)]) ;
% y[k]=q[k]+z[k]
   y(:,k) = q(:,k) + [z(k)+10+v(1,k); z(k)+20+v(2,k); z(k)+30+v(3,k); z(k)+60+v(4,k)]; %+ [v1(:,k); v2(:,k); v3(:,k); v4(:,k)];
end

%plot(z(1,:),'DisplayName','z1')
hold on
plot(y(1,:),'DisplayName','y1')
plot(y(2,:),'DisplayName','y2')
plot(y(3,:),'DisplayName','y3')
plot(y(4,:),'DisplayName','y4')
plot(z+30,'DisplayName','z1 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 
ylim([0,50])
xlim([0 200])

%% case 3: z=[sin(10t)+10+v1 sin(5t)+10+v2 sin(10t)+20+v3 sin(5t)+20+v4... ]'
clc;clear;close all
dt=0.02; 
    t_f=20;
steps=round(t_f/dt);

L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];
q = zeros(8,steps);
y = zeros(8,steps);
L = 1*L;
L_hat = kron(L,eye(2));

y(:,1)=[10 10 20 20 30 30 60 60]';

z1=sin(10*dt*(1:steps));
z2=sin(5*dt*(1:steps));

r=[z1+10 ; z2+10  ; z1+20 ; z2+20  ; z1+30  ; z2+30 ; z1+60 ; z2+60]; %note that when z is a matrix, [z z]≠[z z]'

%generate noise
Ri=0.01;
R = diag([Ri,Ri,Ri,Ri,Ri,Ri,Ri,Ri]);      
v = chol(R)* randn(8,steps);


for k = 2:steps
% q[k]=q[k-1]+dt*(-L)*(q[k-1]-z[k-1])
    %q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + [u(k-1)+10+v(1,k) ; u(k-1)+20+v(2,k) ;u(k-1)+30+v(3,k) ;z1(k-1)+60+v(4,k)]) ;
     q(:,k) = q(:,k-1)+ dt *(-L_hat)*(q(:,k-1) + ( r(:,k-1) + v(:,k-1) ) ) ;
% y[k]=q[k]+z[k]
     y(:,k) = q(:,k) + ( r(:,k) + v(:,k) ); %+ [v1(:,k); v2(:,k); v3(:,k); v4(:,k)];
end

%plot(z(1,:),'DisplayName','z1')
hold on
plot(y(1,:),'DisplayName','y1')
plot(y(2,:),'DisplayName','y2')
plot(y(3,:),'DisplayName','y3')
plot(y(4,:),'DisplayName','y4')
plot(y(5,:),'DisplayName','y1')
plot(y(7,:),'DisplayName','y2')
plot(y(8,:),'DisplayName','y3')

plot(z1+30,'r--','DisplayName','z1 no measurement noise')
plot(z2+30,'g--','DisplayName','z1 no measurement noise')
%plot(zh(2,:),'DisplayName','z2 no measurement noise')
legend 
grid on
ylim([0,50])
xlim([0 200])








