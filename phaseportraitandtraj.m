close all
clear all




%plotAfterPhase()
plotBeforePhase()

%plotDuringPhase()

%plotBeforePhase()


function plotBeforePhase()
    p = getParameter();
    showSmallPhaseSpace()

    p.gamma = 0.21;
    
    g_star = (p.gamma - p.tau*(1 - p.rho))^2 / (p.tau*(1 - p.rho));
    L_star = (p.gamma - p.tau*(1 - p.rho)) / (p.tau*(1 - p.rho));
    scatter([g_star], [L_star], 15, [0 0 1], '^', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)

    


end


function plotDuringPhase()
    p = getParameter();
    showSmallPhaseSpace()
    
    plot(p.a^2*p.tau*(1-p.rho)*[1 1], [p.a/p.theta 20], "LineWidth", 8, "Color", [0 0.6 1])
    scatter([p.a^2*p.tau*(1-p.rho)], [p.a/p.theta], 1e3, [0 0 1], '^', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)



    
end

function plotAfterPhase()
    showSmallPhaseSpace()
    
   
end


% ======== PHASE SPACE ========



function showSmallPhaseSpace()
    fig = figure

    p = getParameter();

    line_width = 1;

    plot([0 3], [p.a/p.theta p.a/p.theta], "LineWidth", line_width, "Color", [0.2 0.8 0.2])
    hold on
    plot(p.tau*p.rho^2/(1-p.rho)*[1 1], [0 20], "LineWidth", line_width, "Color", [0.2 0.2 0.6])
    

    plot([0 p.tau*p.rho^2/(1-p.rho)], p.rho/(1-p.rho)*[1 1], "LineWidth", line_width, "Color", [0.6 0.2 0.2])
    Gs = linspace(p.tau*p.rho^2/(1-p.rho), 3, 100);
    Ls = p.tau*p.rho^2./((1-p.rho)*sqrt(p.tau*(1-p.rho)*Gs));
    plot(Gs, Ls, "LineWidth", line_width, "Color", [0.6 0.2 0.2])

    scatter([p.tau*p.rho^2/(1-p.rho)], [p.a/p.theta], 10, [0.4 1 0.8]/2, 'o', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)
    scatter([p.tau*p.rho^2/(1-p.rho)], [p.rho/(1-p.rho)], 10, [0.8 0.4 0.8]/2, 'o', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 3)


    % solve system
    p.gamma=0.21; % panel b
    p.gamma = p.tau + (p.a*p.tau*(1-p.rho) - p.rho*p.tau); % panel c
    p.gamma = .24; % panel d
    x0 = [1.2632; 768.9279; 11.5427+3]; % G A L
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
    
    col1_A = [255 17 0]/255; % ORIIGNAL
    %col1_A = [255 118 121]/255;
    col1_A = [203 17 0]/255; % 
    col2_A = col1_A;
    
    
    plot(xs(1,:), xs(3,:), '.-', 'MarkerSize', 8, 'LineWidth', 1, 'Color', col1_A)


    box on
    ax = gca;
    ax.FontSize = 6;
    xlim([0 3])
    ylim([0 20])
    width_pt = 1.3*1.3*235/3.12;
    height_pt = 1.3*200*1.5/3;
    
    % Convert points to inches (since MATLAB uses inches for 'Units')
    width_in = width_pt / 72;
    height_in = height_pt / 72;
    
    % Set figure position [left, bottom, width, height]
    set(fig, 'Units', 'inches', 'Position', [1, 1, width_in, height_in]);
    set(gca,'linewidth',1)
end

% ======= HELPERS =======
function p = getParameter()
    p = struct;

    p.rho = 0.879;
    p.tau = 0.15;
    p.a = 11.42;

    p.gamma = p.tau*(1 + p.a)*(1 - p.rho);

    p.c_tilde = 1;
    p.theta = 1;
    p.alpha = 0.6;
end





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
