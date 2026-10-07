% Bifurcation diagram in L
close all
clear all

p = getParameter();
gammas = linspace(0.2, p.gamma, 100);
Ls = (gammas - p.tau*(1-p.rho)) / (p.theta*p.tau*(1-p.rho));

fig = figure; 
plot([0.21 0.21], [7 20], 'Color', 0.4*[1 1 1], 'LineWidth', 1)
hold on
plot(p.gamma*[1 1], [7 20], 'Color', 0.4*[1 1 1], 'LineWidth', 1)
plot(0.24*[1 1], [7 20], 'Color', 0.4*[1 1 1], 'LineWidth', 1)


plot([0.2 0.3], p.a_star*[1 1], 'color', [0.2 0.8 0.2], 'LineWidth', 2, "LineStyle", '--')


plot(gammas, Ls, 'color', "blue", 'LineWidth', 3)
plot([p.gamma p.gamma], [p.a_star 20], 'Color', [0 0.6 1], 'LineWidth', 3)



box on
    ax = gca;
    ax.FontSize = 6;


    xlim([0 3])
    ylim([0 20])


    width_pt = 1.3*1.3*235;
    height_pt = 1.3*200*2/3;
    
    % Convert points to inches (since MATLAB uses inches for 'Units')
    width_in = width_pt / 72;
    height_in = height_pt / 72;
    
    % Set figure position [left, bottom, width, height]
    set(fig, 'Units', 'inches', 'Position', [1, 1, width_in, height_in]);
    
        set(gca,'linewidth',1)
xlim([0.2 0.25])
ylim([7.3 17]);



function p = getParameter()
    p = struct;

    p.rho = 0.879;
    p.tau = 0.15;
    p.a_star = 11.42;

    p.gamma = p.tau*(1 + p.a_star)*(1 - p.rho);

    p.c_tilde = 1;
    p.theta = 1;
    p.alpha = 0.6;
end
