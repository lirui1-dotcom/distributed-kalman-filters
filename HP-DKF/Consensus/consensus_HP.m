

%% Consensus HP Using Ode45
clc; clear; close all;
% initialization ----------------------------------------------------------
% Define the Laplacian matrix L for a square graph of 4 nodes
beta = 0.01;
a=0.00000;

L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];
%L= [2 -1 -1 0;-1 2 0 -1;-1 0 2 -1;0 -1 -1 2];


% L = [1+a -1; -1 1+a];
n = size(L,1);
%L = L+a*eye(n);

fs = 500; % sampling rate

dt=1/fs; % step size
t_f=500;  % simulation time 
steps=round(t_f/dt);

%x0 = [2 3 4 7]'; %initial condition for x
x0 = [1+1 2+1 3+1 6+1]';
%x0 = [0 0 0 0]';
%x0 = [2 6]'; %initial condition for x

x_lim=[0 t_f]; % x axis limit
y_lim=[-5 5]; % y axis limit

% Define the symbolic z(t) functions
syms t_sym
% z_sym = [t_sym, t_sym+0.1, t_sym+1, t_sym+1.1]';
% % f = 10;
% w = 2*pi*f;
w=1;
wL=0.05;
const=1;
%z_sym = [t_sym+1, t_sym+2, t_sym+3, t_sym+6]';
%z_sym = [exp(-const*t_sym)+1, exp(-const*t_sym)+2, exp(-const*t_sym)+3, exp(-const*t_sym)+6]';
z_sym = [cos(w*t_sym)+1, cos(w*t_sym)+2, cos(w*t_sym)+3, cos(w*t_sym)+6]';
%z_sym = [sin(w*t_sym)+sin(wL*t_sym)+2, sin(w*t_sym)+sin(wL*t_sym)+6]';
%z_sym = [sin(w*t_sym)+2, sin(w*t_sym)+6]';
%z_sym = [t_sym, t_sym+2]';
%z_sym = [t_sym+2+sin(0.05*t_sym),t_sym+4+sin(0.05*t_sym)]';



% Take the derivative of z_sym with respect to t_sym
z_dot_sym = diff(z_sym, t_sym);

% Convert symbolic derivative to a function handle
z_dot_fun = matlabFunction(z_dot_sym, 'Vars', t_sym);

% Define the time vector
t = linspace(0, t_f, steps); % simulate from t=0 to tf with n points


% Define the system of equations

x_dot_fun = @(t, x) -beta*L*x - beta*a*eye(n)*x + z_dot_fun(t); % z_dot_fun computes the derivative automatically

% Solve the differential equation using ode45
[t_out, x_out] = ode45(x_dot_fun, t, x0);

% Convert z(t) symbolic functions to function handles for plotting
z_fun = arrayfun(@(z) matlabFunction(z, 'Vars', t_sym), z_sym, 'UniformOutput', false);

% Compute z_avg as the average of z1, z2, z3, and z4
 z_avg = (z_fun{1}(t) + z_fun{2}(t)+ z_fun{3}(t)+ z_fun{4}(t) ) / n;
%z_avg = (z_fun{1}(t) + z_fun{2}(t)) / n;

% Plot the results of x(t) and z_avg on the same plot
figure;
plot(t_out, x_out);%--------------------------------------------------+avg(z0)!!!!!
hold on;
plot(t, z_avg, '--', 'DisplayName', 'z_{avg}(t)');
legend('x_1', 'x_2', 'x_3', 'x_4', 'z_{avg}(t)');
xlabel('Time');
ylabel('x and z_{avg}');
title('Simulation of x\_dot = -Lx + z\_dot with z_{avg}');
xlim(x_lim)
ylim(y_lim)
grid on;

% Plot the original z(t) functions using the generated function handles
figure;
hold on;
for i = 1:length(z_sym)
    plot(t, z_fun{i}(t), 'DisplayName', sprintf('z_%d(t)', i));
end
legend show
xlabel('Time');
ylabel('z(t)');
title('z(t)');
xlim(x_lim)
ylim(y_lim)
grid on;

%figure
%plot(t,z_avg-z_fun{1}(t))
%% Discrete Consensus HP 







%% Frequency response
clear;clc;close all
% Define the Laplacian matrix L for a square graph of 4 nodes
L = [ 2 -1  0 -1;
     -1  2 -1  0;
      0 -1  2 -1;
     -1  0 -1  2];
%L=10*(L+0.01*eye(4));
%L=0.1*(L+50*eye(4));
%L=1*(L+1*eye(4));

%L = 100*[1 -1; -1 1];

% % Number of nodes
% n = 2;
% 
% % Degree of each node
% k = 2;
% 
% % Create a regular graph with 100 nodes and degree 6
% % The adjacency matrix of the graph
% A = zeros(n, n);
% 
% % Construct a degree-6 regular graph
% for i = 1:n
%     for j = 1:k/2
%         A(i, mod(i-1+j, n) + 1) = 1;
%         A(i, mod(i-1-j+n, n) + 1) = 1;
%     end
% end
% 
% % Make sure the adjacency matrix is symmetric
% A = max(A, A');
% 
% % Degree matrix (diagonal matrix with degrees on the diagonal)
% D = diag(sum(A, 2));
% 
% % Laplacian matrix
% L = D - A;

[N,~]=size(L);
% Define the identity matrix I
I = eye(size(L));

% Define the Laplace variable s
s = tf('s');

% Define the transfer function H(s) = s * (sI + L)^-1
H = s /(s*I + L);

% Plot the Bode plot of H(s)
figure;
sigma(H);
grid on;
%set(gca, 'XScale', 'linear');

%xlim([0.1 100]);