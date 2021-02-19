clear
close all

syms k c

% Field

A = 1;
alpha = 0.3;
delta = 0.1;
rho = 0.04;
theta = 0.5;
initial = [1, 3];

tic
%% Forward shooting

% solve for steady state of c and k

REF = [A, alpha, delta, rho, theta];

c_dot = (alpha * k^(alpha - 1) - rho - delta) / theta == 0;
ctemp = @(k) (k^alpha - delta * k);

solutionK = eval(vpasolve(c_dot, k));
solutionC = ctemp(solutionK);

k_0 = solutionK / 2;
c_0 = rand(1);
t = linspace(0, 20, 10000);

% ODE settings

ode_initial1 = [c_0; k_0];
ode_initial2 = [1.5; 10];

option = odeset('NonNegative', [1, 2]);
rck_f = @(t, y) ode(t, y, REF);

[~, Y1] = ode45(rck_f, t, ode_initial1, option);
k_path1 = Y1(:, 2);
c_path1 = Y1(:, 1);

[~, Y2] = ode45(rck_f, t, ode_initial2, option);
k_path2 = Y2(:, 2);
c_path2 = Y2(:, 1);

c_star_k0 = @(k) k.^alpha - delta * k;

% disp(solutionC);
% disp(solutionK);

% Grid search with tolerance, part 1

combine = [c_path1, k_path1];
size = 10000;
grid = linspace(0.63, c_star_k0(k_0), size);
tolerance = 0.01;

for i = grid
    % disp('i:')
    % disp(i);
    true = 0;
    grid_initial = [i; k_0];
    [~, Z] = ode45(rck_f, t, grid_initial, option);
    k_path_grid = Z(:, 2);
    c_path_grid = Z(:, 1);
    if abs(k_path_grid(end) - solutionK) + abs(...
            c_path_grid(end) - solutionC) < tolerance
        true = 1;
    end
    if true == 1
        % disp('The index is:')
        % disp(i);
        break
    end
end

k_0_new = solutionK * 1.5;

% Grid search with tolerance, part 2

grid = linspace(1.5, 1.8, size);
tolerance = 0.01;

for i = grid
    % disp('i:')
    % disp(i);
    true = 0;
    grid_initial=[i; k_0_new];
    [~, ZZ] = ode45(rck_f, t, grid_initial, option);
    k_path_grid_2 = ZZ(:, 2);
    c_path_grid_2 = ZZ(:, 1);
    if abs(k_path_grid_2(end) - solutionK) + abs(...
            c_path_grid_2(end) - solutionC) < tolerance
        true = 1;
    end
    if true == 1
        % disp('The index2 is:')
        % disp(i);
        break
    end
end

% Plot the figure

figure(1)
plot(k_path1, c_path1);
hold on
plot(k_path2, c_path2);
hold on
plot(k_path_grid, c_path_grid);
hold on
plot(k_path_grid_2, c_path_grid_2);
hold on
kt = linspace(0, 20, 10000);

plot(kt, c_star_k0(kt));
hold on
line([solutionK solutionK], [0 10])
axis([0 20 0 2])
title('The Saddle Path (Forward shooting)')
l = legend(...
    'off path sequence 1', 'off path sequence 2', ...
    'saddle path from below', 'saddle path from above', ...
    '$\dot{k}$=0', '$\dot{c}$=0');
set(l, 'interpreter', 'latex')
p = xlabel('capital ($k$)');
set(p, 'interpreter', 'latex')
o = ylabel('consumption ($c$)');
set(o, 'interpreter', 'latex')

%% Backward Integration

% Calculate the Jocobian Matrix and find the eigenvectors

Jocobian = [0 alpha * (alpha - 1) * solutionK^(alpha - 2); -1 rho];
[rightEigenvector, eigenvalueOnDiagonal, leftEigenvector] = eig(Jocobian);

% Calculate the initial point in the direction of eigenvector

epi = 0.01;
below_inital = [solutionC - (rightEigenvector(2, 1) / ...
    (rightEigenvector(2, 1) + rightEigenvector(2, 2))) * 10^(-10); ...
    solutionK - (rightEigenvector(2, 2) / ...
    (rightEigenvector(2, 1) + rightEigenvector(2, 2))) * 10^(-10)];
high_inital = [solutionC + (rightEigenvector(1, 1) / ...
    (rightEigenvector(1, 1) + rightEigenvector(1, 2))) * 10^(-10);
    solutionK + (rightEigenvector(1, 2) / ...
    (rightEigenvector(1, 1) + rightEigenvector(1, 2))) * 10^(-10)];

% Reverse construction

t2 = linspace(100, 0, 8000);

option = odeset('NonNegative', [1, 2]);

[~, BI1] = ode45(rck_f, t2, below_inital, option);
k_path_grid_3 = BI1(:, 2);
c_path_grid_3 = BI1(:, 1);

[~, BI1] = ode45(rck_f, t2, high_inital, option);
k_path_grid_4 = BI1(:, 2);
c_path_grid_4 = BI1(:, 1);

% Plot the figure

figure(2)
plot(k_path_grid_3, c_path_grid_3);
hold on
plot(k_path_grid_4, c_path_grid_4);
hold on
kt = linspace(0, 20, 10000);

plot(kt, c_star_k0(kt));
hold on

line([solutionK solutionK], [0 10])
axis([0 20 0 2])
title('The Saddle Path (Backward Integration)')
l = legend(...
    'saddle path from below', 'saddle path from above', ...
    '$\dot{k}$=0', '$\dot{c}$=0');
set(l, 'interpreter', 'latex')
p = xlabel('capital ($k$)');
set(p, 'interpreter', 'latex')
o = ylabel('consumption ($c$)');
set(o, 'interpreter', 'latex')
toc
