clear all
close all
% =========== PARAMETERS ===========
p = struct;
p.rho   = 0.879;
p.tau   = 0.15;
p.theta = 1;
p.a     = 11.42;
p.alpha = 0.6;
p.gamma = p.tau + (p.a*p.tau*(1-p.rho) - p.rho*p.tau);
p.c_tilde = 1;
% =================================== 8--o
figure(1); hold on


%% FIRST RUN
x0 = [0.048; 0.87; 0.664];
n_iter = 6;
xs = x0;
for i = 1:n_iter
    xs(:, end+1) = F(xs(:, end), p);
end
n_steps = n_iter + 1;
e_series = zeros(1, n_steps);
h_series = zeros(1, n_steps);
z_series = zeros(1, n_steps);
for i = 1:n_steps
    G = xs(1,i); A = xs(2,i); L = xs(3,i);
    e_series(i) = E(G, p);
    h_series(i) = H(G, p);
    z_series(i) = Z(G, A, L, p);
end
t = 0:n_iter;

% save first run for insets later
t1 = t;  xs1 = xs;  e1 = e_series;  h1 = h_series;  z1 = z_series;

col1_A = [0 170 14]/255; % GREEN USED IN OTHER FIGURE
col2_A = col1_A;

subplot(1,3,1)%G
plot(t, xs(1,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
yline(p.tau*p.rho^2/(1-p.rho), 'r--', 'LineWidth', 1.5)
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0 .125])


subplot(1,3,2)%L
plot(t, xs(3,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
yline(p.a/p.theta, 'r--', 'LineWidth', 1.5)
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0 1])


subplot(1,3,3)%A
plot(t, xs(2,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)


ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0.8 1])
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
%x = [0, 45, 45, 0, 0];
%y = [0.75, 0.75, 3, 3, 0.75];
%plot(x, y, 'k-', 'LineWidth', 1);


% Get current figure position ;)
fig = gcf;
pos = fig.Position;

% Hey Shrink height while keeping width and position
new_height = pos(4) * 0.2642;
fig.Position = [pos(1), pos(2) + pos(4) - new_height, pos(3), new_height];

% SAM'S CODE

% ── observables ──────────────────────────────────────────────────────────
function e = E(G, p)
    e = max(0, -p.rho*p.tau + sqrt(p.tau*(1-p.rho)*G));
end

function h = H(G, p)
    e = E(G, p);
    h = (e + p.tau*p.rho) ./ (e + p.tau*p.rho + G);
end

function z = Z(G, A, L, p)
    h = H(G, p);
    z = h.^p.alpha .* (A ./ L).^(1-p.alpha);
end

% ── one step of Lagerlöf's map ────────────────────────────────────────────
function y = F(x, p)
    G = x(1);  A = x(2);  L = x(3);
    e = E(G, p);
    G_next = (e + p.rho*p.tau) * min(p.theta*L, p.a);
    e_next = E(G_next, p);
    A_next = (1 + G_next) * A;
    z = Z(G, A, L, p);
    if z >= p.c_tilde / (1-p.gamma)
        c = p.gamma;
    elseif z >= p.c_tilde
        c = 1 - (p.c_tilde / z);
    else
        c = 0;
        disp("Collapsed case reached")
    end
    L_next = c * L / (p.tau + e_next);
    y = [G_next; A_next; L_next];
end
