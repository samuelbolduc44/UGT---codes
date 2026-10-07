clear all
close all
% =========== PARAMETERS ===========
p = struct;
p.rho   = 0.879;
p.tau   = 0.3;
p.theta = 1;
p.a_star     = 11.42;
p.alpha = 0.6;
p.gamma = 0.245;
p.c_tilde = 1;
p.Am = 10;
% ===================================
figure(1); hold on

%% FIRST RUN
x0 = [0.005; 8.2; 0.025]; % MALTHUSIAN STEADY-STATE
%x0 = [0.005; 8.2; 1.15]; % MODERN-GROWTH STEADY-STATE


n_iter = 200;
xs = x0;
for i = 1:n_iter
    xs(:, end+1) = F(xs(:, end), p);
end
n_steps = n_iter + 1;
e_series = zeros(1, n_steps);
h_series = zeros(1, n_steps);
z_series = zeros(1, n_steps);

e_series = E(xs, p);
h_series = H(xs, p);
z_series = Z(xs, p);

t = 0:n_iter;

% save first run for insets later
t1 = t;  xs1 = xs;  e1 = e_series;  h1 = h_series;  z1 = z_series;

% col1_A = [0 255 17]/255;
col1_A = [0 200 255]/255;
col2_A = [226 134 58]/255; % first one
%col2_A = [100 105 176]/255;
col1_A=col2_A;
subplot(3,2,1)%G
plot(t, xs(1,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
% ylim([0 0.12])
ylim([0 2.5])
ylim([0 0.02])

subplot(3,2,3)%L
plot(t, xs(3,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
% ylim([0 1])
ylim([0 10])
ylim([0 0.25])

subplot(3,2,5)%A
%semilogy(t, xs(2,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
% ONLY USE LOG PLOT FOR THE PURPLE TRAJECTORY
plot(t, xs(2,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
%no ylim assigned for the purple one.
ylim([8.1 10.2])

subplot(3,2,2)
plot(t, e_series, '.-', 'MarkerSize', 8, 'LineWidth', 1,'Color', col2_A);
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
% ylim([0 0.08])
ylim([0 0.1]) %


subplot(3,2,4)
plot(t, h_series, '.-', 'MarkerSize', 8, 'LineWidth', 1,'Color', col2_A)
hold on
ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
ylim([0 1.1])

subplot(3,2,6)
semilogy(t, z_series, '.-', 'MarkerSize', 8, 'LineWidth', 1,'Color', col2_A)
% plot(t, z_series, '.-', 'MarkerSize', 8, 'LineWidth', 1,'Color', col2_A)
hold on

ax = gca; ax.FontSize = 6;
set(gca,'linewidth',1)
% ylim([0.7 10^5])
% ylim([0 3.8])



% ── observables ──────────────────────────────────────────────────────────
function y = F(x, p)
    y = zeros(size(x));

    y(1, :) = Gt(x, p);

    y(2, :) = (1 + y(1, :)).*x(2, :);

    z = Z(x, p);
    if z >= p.c_tilde
        y(3, :) = 4*p.gamma^2/(1+p.gamma)^2 * ((1 - p.c_tilde/z)/(p.tau + E(y, p)).^2).*x(3, :);
    end

end

function g = Gt(x, p)
    % g = (1- exp(-x(3).*(E(x, p) + p.rho*p.tau))) .* (x(2) - p.Am).^2 / (1 + (x(2) - p.Am).^2);
    g = p.a_star*(E(x, p) + p.rho*p.tau)*(1- exp(-x(3)/p.a_star)) .* (x(2) - p.Am).^2 / (1 + (x(2) - p.Am).^2);
    % g = (E(x, p) + p.tau*p.rho).*(x(3) - p.Am).^2 / (1 + (x(3) - p.Am).^2);
end

function e = E(x, p)
    A = (x(1,:) + 2*p.rho*p.tau)/2;
    e = 0.5*(-A + sqrt(A.^2 + (2-4*p.rho/3)*p.tau.*x(1, :) + (4/3)*p.rho*p.tau^2*(1-p.rho)));
end


function z = Z(x, p)
    z = H(x, p).^p.alpha .* (x(2, :) ./ x(3, :)).^(1-p.alpha);
end

function h = H(x, p)
    h = (E(x, p) + p.tau*p.rho/3) ./ (E(x, p) + p.tau*p.rho + x(1, :));
end
