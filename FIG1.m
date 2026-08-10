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
% ===================================
figure(1); hold on

%% FIRST RUN
x0 = [0.048; 0.87; 0.364];
n_iter = 60;
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

col1_A = [255 17 0]/255;
col2_A = [255 127 127]/255;
col2_A=col1_A;

subplot(3,2,1)%G
plot(t, xs(1,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
yline(p.tau*p.rho^2/(1-p.rho), 'r--', 'LineWidth', 1.5)
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0 2.5])
x = [0, 4, 4, 0, 0];
y = [0., 0., .2, .2, 0.];
plot(x, y, 'k-', 'LineWidth', 1);

subplot(3,2,3)%L
plot(t, xs(3,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
yline(p.a/p.theta, 'r--', 'LineWidth', 1.5)
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0 16])
x = [0, 4, 4, 0, 0];
y = [0., 0., 1, 1, 0.];
plot(x, y, 'k-', 'LineWidth', 1);

subplot(3,2,5)%A
semilogy(t, xs(2,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)

subplot(3,2,2)
plot(t, e_series, '.-', 'MarkerSize', 8, 'LineWidth', 1,'Color', col2_A);
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0 0.08])

subplot(3,2,4)
plot(t, h_series, '.-', 'MarkerSize', 8, 'LineWidth', 1,'Color', col2_A)
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0 1.1])

subplot(3,2,6)
semilogy(t, z_series, '.-', 'MarkerSize', 8, 'LineWidth', 1,'Color', col2_A)
hold on
yline(p.c_tilde, 'r--', 'LineWidth', 1.5, 'DisplayName', '\tilde{c}')
yline(p.c_tilde/(1-p.gamma), '--', 'LineWidth', 1.5, ...
      'Color', [1 0.6 0], 'DisplayName', '\tilde{c}/(1-\gamma)')
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0.7 10^5])
x = [0, 45, 45, 0, 0];
y = [0.75, 0.75, 3, 3, 0.75];
plot(x, y, 'k-', 'LineWidth', 1);


%% ── INSETS: zoom into the boxed regions ─────────────────────────────────

% --- Inset for subplot(3,2,1), G, box: x in [0,4], y in [0,0.2] ---
ax_main = subplot(3,2,1);
axPos = get(ax_main,'Position');
insetPos = [axPos(1)+axPos(3)*0.45, axPos(2)+axPos(4)*0.45, ...
            axPos(3)*0.5, axPos(4)*0.5];
axInset1 = axes('Position', insetPos);
plot(axInset1, t1, xs1(1,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold(axInset1,'on')
xlim(axInset1, [0 4])
ylim(axInset1, [0 0.125])
box(axInset1,'on')
set(axInset1,'FontSize',5,'LineWidth',1)

% --- Inset for subplot(3,2,3), L, box: x in [0,4], y in [0,1] ---
ax_main = subplot(3,2,3);
axPos = get(ax_main,'Position');
insetPos = [axPos(1)+axPos(3)*0.45, axPos(2)+axPos(4)*0.45, ...
            axPos(3)*0.5, axPos(4)*0.5];
axInset2 = axes('Position', insetPos);
plot(axInset2, t1, xs1(3,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold(axInset2,'on')
xlim(axInset2, [0 4])
ylim(axInset2, [0 .75])
box(axInset2,'on')
set(axInset2,'FontSize',5,'LineWidth',1)

% --- Inset for subplot(3,2,6), z, box: x in [0,45], y in [0.75,3] ---
ax_main = subplot(3,2,6);
axPos = get(ax_main,'Position');
insetPos = [axPos(1)+axPos(3)*0.65, axPos(2)+axPos(4)*0.65, ...
            axPos(3)*0.65, axPos(4)*0.65];
axInset3 = axes('Position', insetPos);
semilogy(axInset3, t1, z1, '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col2_A)
hold(axInset3,'on')
xlim(axInset3, [0 45])
ylim(axInset3, [0.75 3])
box(axInset3,'on')
set(axInset3,'FontSize',5,'LineWidth',1)
yline(p.c_tilde, 'r--', 'LineWidth', 1.5, 'DisplayName', '\tilde{c}')
yline(p.c_tilde/(1-p.gamma), '--', 'LineWidth', 1.5, ...
      'Color', [1 0.6 0], 'DisplayName', '\tilde{c}/(1-\gamma)')

% Use vector graphics renderer
set(gcf,'renderer','painters')

axs = findall(gcf,'Type','Axes');
for k = 1:numel(axs)
    axs(k).DataAspectRatioMode = 'auto';
    axs(k).PlotBoxAspectRatio = [1.25*1.61803398875 1 1];
end
for k = 1:numel(axs)
    pos = axs(k).Position;
    pos(3) = pos(3)*1.10;
    pos(4) = pos(4)*1.15;
    pos(1) = pos(1)-0.01;
    pos(2) = pos(2)-0.01;
    axs(k).Position = pos;
end
set(gcf,'Units','points', 'Position',[100,100,625,425])
pos = get(gcf, 'Position');
set(gcf, 'Units', 'points', 'Position', [pos(1), pos(2), pos(3)*0.75, pos(4)*0.75])

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