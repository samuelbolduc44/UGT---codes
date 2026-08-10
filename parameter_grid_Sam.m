clear all
close all

%% parameters
p = struct;
p.rho     = 0.879;
p.tau     = 0.15;
p.theta   = 1;
p.a       = 11.42;
p.alpha   = 0.6;
p.gamma   = p.tau + (p.a*p.tau*(1-p.rho) - p.rho*p.tau);
p.c_tilde = 1;

x0 = [0.048; 0.87; 0.364]; % G0 A0 L0   (e0 removed: e is determined by G)
%x0 = [0.075; 0.87; 0.364]; % G0 A0 L0
%x0 = [0.125; 0.87; 0.364]; % G0 A0 L0

%% Grid definition
n_A = 10;
n_L = 10;
A_grid = linspace(0, 25, n_A);
L_grid = linspace(0., 25, n_L);

%A_grid = linspace(0, 5, n_A);
%L_grid = linspace(0, 5, n_L);

%A_grid = linspace(3., 5, n_A);
%L_grid = linspace(0., .2, n_L);

%% Parallel sweep
collapse_matrix = NaN(n_A, n_L); %set matrix

parfor i = 1:n_A
    col = NaN(1, n_L);
    for j = 1:n_L
        x0_ij = x0;
        x0_ij(2) = A_grid(i); %update A0
        x0_ij(3) = L_grid(j); %.....  L0
        [~, col(j)] = run(p, x0_ij);
    end
    collapse_matrix(i, :) = col;
end

%% Display result
%disp(collapse_matrix)

%% Visualize
fig = figure(1)
hold on

% Colormap from https://colorhunt.co/palette/213c516594b1ddaed3eeeeee
colorNaN = [238 238 238]/255;

% Base colors
baseColors = [101/255 148/255 177/255;     % itr = 1
              221/255 174/255 211/255];    % itr = 2

% MATLAB default colors for cycling (7 colors)
matlabColors = [0    0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.9290 0.6940 0.1250;
                0.4940 0.1840 0.5560;
                0.4660 0.6740 0.1880;
                0.3010 0.7450 0.9330;
                0.6350 0.0780 0.1840];

% Desired number of iterations
nIter = 60; % change this to whatever you need

% Create coloritr for any nIter
coloritr = [baseColors; 
            repmat(matlabColors, ceil((nIter-2)/size(matlabColors,1)), 1)];
coloritr = coloritr(1:nIter, :); % Trim to exact size needed

%%
% Marker size in data units
dA = A_grid(2) - A_grid(1);
dL = L_grid(2) - L_grid(1);

for i = 1:n_A
    for j = 1:n_L
        val = collapse_matrix(i, j);
        if isnan(val)
            c = colorNaN;
        elseif val <= size(coloritr, 1)
            c = coloritr(val, :);
        else
            c = coloritr(end, :); % clamp to last color if itr > 4
        end
        rectangle('Position', [L_grid(j)-dL/2, A_grid(i)-dA/2, dL, dA], ...
                  'FaceColor', c, 'EdgeColor', 'none')
    end
end

%xlabel('$L_0$', 'Interpreter', 'latex')
%ylabel('$A_0$', 'Interpreter', 'latex')
%title('Projection phase space analysis: e0=1, kappa=11, beta=1')
%ax = gca; ax.FontSize = 14; ax.TickLabelInterpreter = 'latex';
xlim([L_grid(1)-dL/2, L_grid(end)+dL/2])
ylim([A_grid(1)-dA/2, A_grid(end)+dA/2])

width_pt = 1.3*1.3*235;
height_pt = 1.3*200;

% Convert points to inches (since MATLAB uses inches for 'Units')
width_in = width_pt / 72;
height_in = height_pt / 72;

% Set figure position [left, bottom, width, height]
set(fig, 'Units', 'inches', 'Position', [1, 1, width_in, height_in]);


ax = gca;
set(gca,'linewidth',1)
ax.XTickLabel=[];
ax.YTickLabel=[];
%print(fig, 'rectangle.png', '-dpng', '-r200');

%% Functions

function [collapsed, collapse_iter] = run(p, x0)
    n_iter = 60;
    xs = x0;
    collapsed = false;
    collapse_iter = NaN;
    %xprev = x0;

    for i = 1:n_iter
        %x_prev = xs(:,end);
        [x_next, col] = F(xs(:, end), p);
        xs(:, end+1) = x_next;
        if col
            collapsed = true;
            collapse_iter = i;
            break
        end
    end
end

function [y, collapsed] = F(x, p)
    G = x(1);  A = x(2);  L = x(3);

    e = E(G, p);                                  % observable, from G

    G_next = (e + p.rho*p.tau) * min(p.theta*L, p.a);
    e_next = E(G_next, p);                        % needed for L_next
    A_next = (1 + G_next) * A;

    z = Z(G, A, L, p);
    collapsed = false;
    if z > p.c_tilde / (1 - p.gamma)
        c = p.gamma;
    elseif z > p.c_tilde
        c = 1 - (p.c_tilde / z);
    else
        c = 0;
        collapsed = true;
    end

    L_next = c * L / (p.tau + e_next);

    y = [G_next; A_next; L_next];
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

function color = getIterationColor(itr)
    % Fixed colors
    colorBase = [101 148 177;     % itr = 1
                 221 174 211]/255; % itr = 2
    
    % Basic MATLAB colors
    matlabColors = [0    0.4470 0.7410;
                    0.8500 0.3250 0.0980;
                    0.9290 0.6940 0.1250;
                    0.4940 0.1840 0.5560;
                    0.4660 0.6740 0.1880;
                    0.3010 0.7450 0.9330;
                    0.6350 0.0780 0.1840];
    
    if itr == 1
        color = colorBase(1,:);
    elseif itr == 2
        color = colorBase(2,:);
    else
        % Cycle through matlabColors for itr >= 3
        idx = mod(itr-3, size(matlabColors,1)) + 1;
        color = matlabColors(idx, :);
    end
end