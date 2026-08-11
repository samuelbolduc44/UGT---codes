



% plotBeforePhase()

% plotDuringPhase()

plotAfterPhase()


function plotBeforePhase()
    p = getParameter();
    showSmallPhaseSpace()

    p.gamma = 0.21;
    
    g_star = (p.gamma - p.tau*(1 - p.rho))^2 / (p.tau*(1 - p.rho));
    L_star = (p.gamma - p.tau*(1 - p.rho)) / (p.tau*(1 - p.rho));
    scatter([g_star], [L_star], 1e3, [0 0 1], '^', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)

    exportgraphics(gcf,'before_phase.eps','ContentType', 'vector')
    close

end


function plotDuringPhase()
    p = getParameter();
    showSmallPhaseSpace()
    
    plot(p.a_star^2*p.tau*(1-p.rho)*[1 1], [p.a_star/p.theta 20], "LineWidth", 8, "Color", [0 0.6 1])
    scatter([p.a_star^2*p.tau*(1-p.rho)], [p.a_star/p.theta], 1e3, [0 0 1], '^', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)

    exportgraphics(gcf,'during_phase.eps','ContentType', 'vector')
    close

    
end

function plotAfterPhase()
    showSmallPhaseSpace()

    % exportgraphics(gcf,'after_phase.eps','ContentType', 'vector')
    % close
end


% ======== PHASE SPACE ========

function showPhaseSpace()
    figure

    p = getParameter();

    plot([0 3], [p.a_star/p.theta p.a_star/p.theta], "LineWidth", 2, "Color", [0.2 0.8 0.2])
    hold on
    plot(p.tau*p.rho^2/(1-p.rho)*[1 1], [0 20], "LineWidth", 2, "Color", [0.2 0.2 0.6])
    

    plot([0 p.tau*p.rho^2/(1-p.rho)], p.rho/(1-p.rho)*[1 1], "LineWidth", 2, "Color", [0.6 0.2 0.2])
    Gs = linspace(p.tau*p.rho^2/(1-p.rho), 3, 100);
    Ls = p.tau*p.rho^2./((1-p.rho)*sqrt(p.tau*(1-p.rho)*Gs));
    plot(Gs, Ls, "LineWidth", 3, "Color", [0.6 0.2 0.2])

    scatter([p.tau*p.rho^2/(1-p.rho)], [p.a_star/p.theta], 1e2, [0.4 1 0.8]/2, 'o', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)
    scatter([p.tau*p.rho^2/(1-p.rho)], [p.rho/(1-p.rho)], 1e2, [0.8 0.4 0.8]/2, 'o', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)

    box on
    ax = gca;
    set(gca,'linewidth', 1.5)

    x0 = 500;
    y0 = 200;
    width = 1200;
    height = 1200;
    set(gcf,'position',[x0,y0,width,height])

    xticks([0 1 2 3])
    yticks([0 10 20])

    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);

    xlim([0 3])
    ylim([0 20])
    pbaspect([1 1 1])

    

end


function showSmallPhaseSpace()
    figure

    p = getParameter();

    line_width = 10;

    plot([0 3], [p.a_star/p.theta p.a_star/p.theta], "LineWidth", line_width, "Color", [0.2 0.8 0.2])
    hold on
    plot(p.tau*p.rho^2/(1-p.rho)*[1 1], [0 20], "LineWidth", line_width, "Color", [0.2 0.2 0.6])
    

    plot([0 p.tau*p.rho^2/(1-p.rho)], p.rho/(1-p.rho)*[1 1], "LineWidth", line_width, "Color", [0.6 0.2 0.2])
    Gs = linspace(p.tau*p.rho^2/(1-p.rho), 3, 100);
    Ls = p.tau*p.rho^2./((1-p.rho)*sqrt(p.tau*(1-p.rho)*Gs));
    plot(Gs, Ls, "LineWidth", line_width, "Color", [0.6 0.2 0.2])

    scatter([p.tau*p.rho^2/(1-p.rho)], [p.a_star/p.theta], 5e2, [0.4 1 0.8]/2, 'o', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)
    scatter([p.tau*p.rho^2/(1-p.rho)], [p.rho/(1-p.rho)], 5e2, [0.8 0.4 0.8]/2, 'o', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)

    box on
    ax = gca;
    set(gca,'linewidth', 3)

    x0 = 500;
    y0 = 200;
    width = 1200;
    height = 1200;
    set(gcf,'position',[x0,y0,width,height])

    xticks([0 1 2 3])
    yticks([0 10 20])

    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);

    xlim([0 3])
    ylim([0 20])
    pbaspect([1 1 1])

    

end

% ======= HELPERS =======
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
