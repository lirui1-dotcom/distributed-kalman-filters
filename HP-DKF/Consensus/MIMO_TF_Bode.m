% Define the transfer function
s = tf('s');

% L = [ 2 -1 -1 0 ;
%      -1  2  0 -1;
%      -1  0  2 -1;
%       0 -1 -1  2];
L = [2 -1; -1 2];
% Define the identity matrix I
I = eye(size(L));

% Define the transfer function H(s) = s * (sI + L)^-1
H = s /(s*I + L);

% Frequency range
omega = logspace(-1, 2, 500);

% Number of singular values to plot (based on the size of L)
numSingularValues = size(L, 1);

% Preallocate storage for singular values
singularValues = zeros(numSingularValues, length(omega));

% Compute singular values at each frequency
for i = 1:length(omega)
    H_omega = evalfr(H, 1j*omega(i));
    [~, S, ~] = svd(H_omega);
    singularValues(:, i) = diag(S);
end

% Plot the results
figure;
hold on;

% Define a set of colors and line styles to cycle through
colors = {'b-', 'r--', 'g-.', 'm:', 'c-', 'k--', 'y-.', 'b:'};  

% Plot each singular value
for j = 1:numSingularValues
    semilogx(omega, 20*log10(singularValues(j, :)), colors{mod(j-1, length(colors)) + 1}, 'LineWidth', 1.5);
end

grid on;
xlabel('Frequency (rad/s)');
ylabel('Singular Values (dB)');
title('Singular Value Decomposition of the Transfer Function');
legend(arrayfun(@(x) sprintf('Singular Value %d', x), 1:numSingularValues, 'UniformOutput', false));
hold off;