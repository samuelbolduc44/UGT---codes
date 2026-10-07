clear all
close all
% 0    0.4470 0.7410
Ls = [0 2];
As = [0 15];


G0 = 0.005;

% G0 = 0.0;

tic
points = phaseSweep(Ls, As, G0);
fig=plotSweep(points, Ls, As, G0 )
toc

% debug()

width_pt = 1.3*1.3*235;
height_pt = 1.3*200;

% Convert points to inches (since MATLAB uses inches for 'Units')
width_in = width_pt / 72;
height_in = height_pt / 72;

% Set figure position [left, bottom, width, height]
set(fig, 'Units', 'inches', 'Position', [1, 1, width_in, height_in]);
%set(gcf,'renderer','painters')

ax = gca;
set(gca,'linewidth',1)
ax = gca; ax.FontSize = 6;
ax.XTickLabel=[];
ax.YTickLabel=[];



print(fig, 'rectangle.png', '-dpng', '-r200')

% ======== PHASE SWEEP FUNCTIONS ========
function fig = plotSweep(points, xlimits, ylimits, G0)

    grid_points = getGrid(xlimits, ylimits, G0, 400);
    X = grid_points(:, :, 3);
    Y = grid_points(:, :, 2);
    ps = [X(:)'; Y(:)'];

    ps(:, 32)
    points(:, 32)

    % process colors
    colours = zeros(1, size(ps, 2));

    colours(points(2, :) < 11) = 2;
    colours(points(3, :) == 0) = 1;
    
    colours = colours + 1;

    c_map = 1/255*[100 176 109;
            101 148 177; 
            176 69 78];

    colours(32)
    colours(1) = 1;
    % clim([1 3])
    

    fig = figure 
    scatter(ps(1, :), ps(2, :), 10, colours, 'filled', 'square')
    hold on

    p = getParam();
    h = H([G0; 0; 0], p);
    grad = h^(-p.alpha/(1-p.alpha));
    %plot([0 5], [0, grad*5], 'LineWidth', 2)
    
    colormap(c_map)
    % circle(0,0,1);

    box on
    set(gca,'linewidth',1.5)
    %plot(0.025,8.2,'*')
    %plot(1.15,8.2,'*')
    xlim(xlimits)
    ylim(ylimits)

    
end



function debug()
    X0 = [0.048; 3.88888; 0.555];

    p = getParam();

    Z(X0, p)
    F(X0, p)

end


% ======== PHASE SWEEP FUNCTIONS ========

function points = phaseSweep(xlimits, ylimits, G0)
    p = getParam();

    n = 400; % 2000
    grid = getGrid(xlimits, ylimits, G0, n);
    Gs = grid(:, :, 1);
    As = grid(:, :, 2);
    Ls = grid(:, :, 3);

    points = [Gs(:)'; As(:)'; Ls(:)'];
    
    
    for i = 1:300
        points = F(points, p);
    end
    % points(:, 18)

end




function points = getGrid(xlimits, ylimits, G0, n)
    % [R, Theta] = meshgrid(linspace(0, 1, n), linspace(0, 2*pi, n));
    % X = R .* cos(Theta);
    % Y = R .* sin(Theta);

    [X, Y] = meshgrid(linspace(xlimits(1), xlimits(2), n), linspace(ylimits(1), ylimits(2), n));
    % [X, Y] = meshgrid(linspace(0.2, 3, n), linspace(-1.6, -1.1, n));

    
    points(:, :, 3) = X;
    points(:, :, 2) = Y;
    points(:, :, 1) = G0 * ones(size(points(:, :, 2)));

end



% ======== SYSTEM FUNCTIONS ========

function y = F(x, p)
    y = zeros(size(x));

    y(1, :) = G(x, p);

    y(2, :) = (1 + y(1, :)).*x(2, :);

    z = Z(x, p);
    y(3, z >= p.c_tilde) = 4*p.gamma^2/(1+p.gamma)^2 * ((1 - p.c_tilde./z(z >= p.c_tilde))./(p.tau + E(y(:, z >= p.c_tilde), p))).^2.*x(3, z >= p.c_tilde);

end


function g = G(x, p)
    g = p.a_star*(E(x, p) + p.rho*p.tau).*(1- exp(-x(3, :)/p.a_star)) .* (x(2, :) - p.Am).^2 ./ (1 + (x(2, :) - p.Am).^2);
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


function p = getParam()
    p = struct;
    p.rho = 0.879;
    p.tau = 0.3;
    p.theta = 1;
    p.gamma = 0.245;
    p.c_tilde = 1;
    p.alpha = 0.6;
    p.a_star = 11.42;
    p.X = 1;
    p.Am = 10;

    % p = struct;
    % p.rho = 0.6545; % 0.63
    % p.tau = 0.15;
    % p.theta = 1;
    % p.gamma = 0.126761; % tau + (a*tau*(1-rho) - rho*tau); 
    % p.c_tilde = 1;
    % p.alpha = 0.6;
    % p.a_star = 50;
    % p.X = 1;
    % 
    % p.Am = 10;
end

