% Bifurcation diagram in L

p = getParameter();
gammas = linspace(0.2, p.gamma, 100);
Ls = (gammas - p.tau*(1-p.rho)) / (p.theta*p.tau*(1-p.rho));


plot([0.21 0.21], [7 20], 'Color', 0.4*[1 1 1], 'LineWidth', 1)
hold on
plot(p.gamma*[1 1], [7 20], 'Color', 0.4*[1 1 1], 'LineWidth', 1)
plot(0.24*[1 1], [7 20], 'Color', 0.4*[1 1 1], 'LineWidth', 1)


plot([0.2 0.3], p.a_star*[1 1], 'color', [0.2 0.8 0.2], 'LineWidth', 2, "LineStyle", '--')


plot(gammas, Ls, 'color', "blue", 'LineWidth', 3)
plot([p.gamma p.gamma], [p.a_star 20], 'Color', [0 0.6 1], 'LineWidth', 3)





box on
ax = gca;
set(gca,'linewidth', 1.5)

x0 = 500;
y0 = 200;
width = 1200;
height = 1200;
set(gcf,'position',[x0,y0,width,height])

xticks([0.2 0.225 0.25])
yticks([10 15])

set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

xlim([0.2 0.25])
ylim([7.3 17]);
pbaspect([3 1 1])

exportgraphics(gcf,'bifurcation.eps','ContentType', 'vector')
close

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