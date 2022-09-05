%% Multi-Marginal Optimal mass transport problem
% This script numerically solves an OMT problem as follows:
% 1. Transportation cost: C_{ijk} = \frac{x_{i}^2}{t_{k}} + \frac{y_{i}^2}{T - t_{k}}
% 2. Transportation map: Pi_{ijk} with its three marginials: p_i, q_i, r_k
%    such that \sum_{jk} = p_{i}, \sum_{ik} = q_{j}, \sum_{ij} = r_k.
%    p_i and q_j are two probability vector   
% 3. The mariginal along the t cooridinates has an upper-bound cap*T/n holds all the time.
% 4. p is transported to q in the time interval [0,T]

% The optimizer is the vector r_k which minimize the transportation cost:
% min_{ijk} C_{ijk}.*\Pi_{ijk}

% The code works in MATLAB 2019 a,b
% CVX with MOSEK solver is used:
% Installation of CVX: http://cvxr.com/cvx/

%% This Script is for Dirac transportation and plot, check note for the desired figure we want to see 
clear all
clc
%% Define the two probability distributions p_{i} and q_{j} and cooridinates X(p), Y(q)
% Consider p_{i} and q_{j} are two discretes distributions 
% p_{i}: i = 1, 2,..., n1
% q_{j}: j = 1, 2,..., n2
 
%% Distributions generator in 2-D
% Dimensions n0,n1 of the clouds
n0 = 40;
n1 = 40;


gauss = @(q,a,c)a*randn(2,q)+repmat(c(:), [1 q]);

% Let X0 between [-1,0]
X0 = randn(2,n0)*.3;
X0(1,:) = rescale(X0(1,:),-1,-0.1);
X0(2,:) = rescale(X0(2,:),-1,-0.1);

x = -sqrt(X0(1,:).^2+X0(2,:).^2);

% Let X1 between [0,1]
X1 = [gauss(n1/2,.5, [0 1.6]) gauss(n1/4,.3, [-1 -1]) gauss(n1/4,.3, [1 -1])];
X1(1,:) = rescale(X1(1,:),0.1,0.9);
X1(2,:) = rescale(X1(2,:),0.1,0.9);
y = -sqrt(X1(1,:).^2+X1(2,:).^2);

normalize = @(a)a/sum(a(:));

p0 = normalize(rand(n0,1));
p1 = normalize(rand(n1,1));

figure(1)
figure(1)
subplot(2,1,1);
bar(x, p0, 'k'); axis tight;
title('\rho_{0}')
subplot(2,1,2);
bar(y, p1, 'k'); axis tight;
title('\rho_{1}')

%%
myplot = @(x,y,ms,col)plot(x,y, 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 1);

figure(2)
hold on;
% two dimension case
for i=1:length(p0)
    myplot(X0(1,i), X0(2,i), p0(i)*length(p0)*5, 'b');
end

for i=1:length(p1)
    myplot(X1(1,i), X1(2,i), p1(i)*length(p1)*5, 'r');
end

% diagonal mass summation
for i=1:length(p0)
    myplot(X0(1,i), -X0(1,i), p0(i)*length(p0)*5, 'b');
end

for i=1:length(p0)
    myplot(X1(1,i), -X1(1,i), p1(i)*length(p1)*5, 'r');
end

% Original point, location of the 'toll'
% scatter(0,0,,'filled')

% % rho-0 correspondence
% plot(X0(1,:),X0(1,:),X0(2,:),-X0(1,i),'b')
% plot(X0(1,:),-X0(1,:),'b','LineWidth',1);
% % rho-1 correspondence
% plot(X1(1,:),-X1(1,:),'r','LineWidth',1);

% Projection rho 0
for i=1:length(p0)
    line([X0(1,i),X0(1,i)],[X0(2,i),-X0(1,i)],'Color','k');
end

% Projection rho 1
% Projection rho 0
for i=1:length(p1)
    line([X1(1,i),X1(1,i)],[X1(2,i),-X1(1,i)],'Color','k');
end


axis([-1 1 -1 1]); 
axis off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off

title('2-D')
%% Define the time t_{i} \in [0, T] (tk):
% Discretization of the time mariginal by nt (not exceed 500)
nt = 200;

% Define the time interval [0,T], notice that T is not neccesary to be 1.
T = 1;
t = T*(0+0.01:nt-0.01)'/nt;

%% Define the cost function C(x,y,t):
% C(x,y,t) = C_{ijk} = (x_i)^2./t_k + (y_i)^2/(T - t_k)
% Compute the 3-dimensional Cost matrix
[X,Y,tk] = meshgrid(x,y,t);
C = X.^2./tk + Y.^2./(T-tk);

% The order of the x-y axis need to be changed by the permute function,
Cost = permute(C,[2,1,3]);

%% CVX solver 
% The optimization problem:
% Variables: 
% The multi marginial coupling \Pi, the 'toll' constraint on time r_{k}
% Objective fucntion: \sum_{i}\sum_{j}\sum_{k} Cost(i,j,k).*\Pi 
% Constraints: 
% 1) marginal p_i, q_j, r_k
% 2) lower bound of \Pi_{i,j,k}: \Pi_{i,j,k} \geq 0
% 3) lower and upper bound of r_{k}: 0 \geq r_{k} \leq (capacity.*T./nt)

% The CVX solver: MOSEK needs to be installed. 
% The sum function in the CVX package do not support high dimension summation
% The capacity constraint should be larger than 1: cap*T/nt>1

% Define the capacity 
cap = 2.5;
capacity = cap*T/nt;
%%
cvx_begin
cvx_solver mosek
cvx_precision best
    variable Pi(n1,n2,nt)
    variable rk(nt)
    minimize(sum(sum(sum(Cost.*Pi))))
    subject to  
    % Define multi-marginal coupling counstraints:
    % 1.Marginal constraints on x coordinate
      squeeze(sum(sum(permute(Pi,[2,3,1])))) == p;
    % 2.Marginal constraints on y coordinate   
      squeeze(sum(sum(permute(Pi,[3,1,2])))) == q;
    % 3.Marinigal on t-coordinate
      squeeze(sum(sum(Pi))) == rk;  
    % 4.Element-wise lower bound for \Pi
      Pi >= 0;
    % 5.Capacity/Toll constraint
      0 <= rk <= capacity;    % r_{k} with upperbound,
      % 0 <= rk               % r_{k} without upperbound.
cvx_end


%% Check 3 marginal constraints:
fprintf('Constraints deviation (should be 0):'); 
sum(squeeze(sum(sum(permute(Pi,[2,3,1])))) - p)
sum(squeeze(sum(sum(permute(Pi,[3,1,2])))) - q)
sum(squeeze(sum(sum(Pi)))- rk)

%% Show that t marginal is a prbability measure
fprintf('Time marginal measure check:'); 
sum(rk)

%% Plots of rk, evolutionary probability distributions at moment t
figure(3)
flux = plot(t,rk);
flux.LineWidth = 2;
flux.Color ='#0072BD';
title('Marginal of Time between [0,T]')

figure(4)
bar(t,rk, 'k'); axis tight;
title('Marginal of Time between [0,T]')

%% 3D plot with 3 marginals
% Things need to be plotted here:
% 1. p  with x-axis, 2. q with y-axis, 3. rk with t-axis, 4. capacity 
% Original coordinates:
% x = (0:n1-1)'/n1-1;         % x \in [-1,0]
% y = (1:n2)'/n2;             % y \in [0,1]
% t = T*(0+0.01:nt-0.01)'/nt; % t \in [0,1]

% Generate new coordinates for 3D plot
xx1 = x;
xx2 = zeros(size(xx1));
xx = [xx1,xx2,p];
yy1 = y;
yy2 = ones(size(yy1));
yy = [yy1,yy2,q];
tt1 = t;
tt2 = zeros(size(tt1));
tollheight = capacity*ones(size(tt1));


%% figure: Line plot3
figure(5)
p1 = plot3(xx1,xx2,p,'LineWidth',3);
hold on
p2 = plot3(yy1,yy2,q,'LineWidth',3);
scaler = 20;
% tplot = plot(rk*scaler,tt1,'LineWidth',2,'Color','k','LineStyle',':');
tplot = plot3(tt2,tt1,rk,'LineWidth',2,'Color','k');
toll = plot3(tt2,tt1,tollheight,'LineWidth',4,'Color','r','LineStyle',':');
hold off
grid on
xlabel('X/Y')
ylabel('Time')
zlabel('Probability Density')
title('3-D Line Plot')


%% Coupling of p and q illustration
pqCoupling = sum(Pi,3);
ptCoupling = squeeze(sum(Pi,2));
qtCoupling = squeeze(sum(Pi,1));

figure(6)
plot(t,rk)
%% Display the optimal assignement

figure(7)
p1 = plot3(xx1,xx2,p,'LineWidth',1.5,'Color','k');

hold on
p2 = plot3(yy1,yy2,q,'LineWidth',1.5,'Color','k');
tplot = plot3(tt2,tt1,rk,'LineWidth',1.5,'Color','k');
toll = plot3(tt2,tt1,tollheight,'LineWidth',1.5,'Color','r','LineStyle',':');
% bar3(tt2,tt1,tollheight);

% Plot the grey 'area' by Line
% pt Coupling (correct):
[Ipt,Jpt,~] = find(ptCoupling);
for ipt=1:length(Ipt)
    pt_assignement = plot([xx(Ipt(ipt),1),tt2(Jpt(ipt))],[xx(Ipt(ipt),2),tt1(Jpt(ipt))],'LineWidth',3,'Color',[0.8 0.8 0.8]);
end

% qt Coupling (order fixed):
[Iqt,Jqt,~] = find(qtCoupling);
for iqt=1:length(Iqt)
    qt_assignement = plot([yy(Iqt(iqt),1),tt2(Jqt(iqt))],[yy(Iqt(iqt),2),tt1(Jqt(iqt))],'LineWidth',3,'Color',[0.8 0.8 0.8]);
end

% Plot the Dash line within the grey area to show point-wise assigments
for iptdash=1:10:length(Ipt)
    pt_assignement_dash = plot([xx(Ipt(iptdash),1),tt2(Jpt(iptdash))],[xx(Ipt(iptdash),2),tt1(Jpt(iptdash))],'LineWidth',0.5,'LineStyle','--','Color','k');
end

for iqtdash=1:10:length(Iqt)
    qt_assignement_dash = plot([yy(Iqt(iqtdash),1),tt2(Jqt(iqtdash))],[yy(Iqt(iqtdash),2),tt1(Jqt(iqtdash))],'LineWidth',0.5,'LineStyle','--','Color','k');
end

% Color the 'bumps' with dark grey
% Bump: 1,2,3,... need manually input y coordinates of the bumps

yarea = [0.095 0.105 0.355 0.365 0.47 0.48 0.51 0.545 0.595 0.615 0.79 0.805];
xarea = zeros(size(yarea));

zperiod = [0 capacity capacity 0];
zarea = [zperiod,zperiod,zperiod];

tollgrey = fill3(xarea,yarea,zarea, [0.2 0.2 0.2]);
tollgrey.FaceAlpha = 0.5;


hold off
xlabel('Space','FontSize',14)
ylabel('Time','FontSize',14)
zlabel('Probability Density','FontSize',14)