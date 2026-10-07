close all
clear all
%% import boundary of M2
%bd = load('0.048.mat'); % a
bd = load('0.35.mat'); % b
%bd = load('0.42.mat'); % c
x_coord = bd.L0;
y1_coord = bd.branch1;
y2_coord = bd.branch2;
fig = figure; hold on
plot(x_coord,y1_coord,'-','LineWidth',1)
plot(x_coord,y2_coord,'-','LineWidth',1)
%% boundary of M1
% parameters
p = struct;
p.rho   = 0.879;
p.tau   = 0.15;
p.theta = 1;
p.a     = 11.42;
p.alpha = 0.6;
p.gamma = p.tau + (p.a*p.tau*(1-p.rho) - p.rho*p.tau);
p.c_tilde = 1;
% variables
syms xvar
g_0 = 0.35; % TO MANUALLY ADJUST
e_0=0;
fplot(xvar*(p.c_tilde*(1+g_0/(e_0+p.tau*p.rho))^p.alpha)^(1/(1-p.alpha)),'b')
%% plotting
axis([0 25 0 250])
width_pt = 1.3*1.3*235*0.2791;
height_pt = 1.3*200*0.4269;

% Convert points to inches (since MATLAB uses inches for 'Units')
width_in = width_pt / 72;
height_in = height_pt / 72;

% Set figure position [left, bottom, width, height]
set(fig, 'Units', 'inches', 'Position', [1, 1, width_in, height_in]);

ax = gca;
set(gca,'linewidth',1)
ax.XTickLabel=[];
ax.YTickLabel=[];
set(gcf,'renderer','painters');
