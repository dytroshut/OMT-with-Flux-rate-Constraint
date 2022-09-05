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
 
%% Random distribution generator
% 
n1 = 200;  % (not exceed 1000)
n2 = 200;  % (not exceed 1000)

x = (0:n1-1)'/n1-1; % x \in [-1,0]
y = (1:n2)'/n2;     % y \in [ 0,1]

% Useful functions 
normalize = @(a)a/sum(a(:));                           % Normalization function 
Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) ); % Gaussian distribution generator

sigma0 = 0.01; sigma1 = .02; sigma2 = .03; sigma3 = .04; sigma4 = .05;

% % Case 1:
% % Marginal p
% p = Gaussian(x, -.6,sigma4) + 3.*Gaussian(x, -.7,sigma2);
% % Marginal q
% q = 0.2.*Gaussian(y, .4,sigma4)+ 4.*Gaussian(y, .4,sigma2)+ 1.5*Gaussian(y, .5,sigma2)+ 1*Gaussian(y, .6,sigma2)+ 5.*Gaussian(y, .7,sigma3);

% Case 2:
% % Marginal p
p = 5.*Gaussian(x, -.77,sigma3) + 1.5.*Gaussian(x, -.7,sigma0) + 2.5.*Gaussian(x, -.3,sigma2) + 4.*Gaussian(x, -.2,sigma3);
% % Marginal q
q = 0.2.*Gaussian(y, .2,sigma4)+ 3.*Gaussian(y, .3,sigma2)+ 1.5*Gaussian(y, .4,sigma2)+ 1*Gaussian(y, .7,sigma2)+ 5.*Gaussian(y,.8,sigma3);

p = normalize(p);
q = normalize(q);

% Plot the two distributions p and q
figure(1)
subplot(2,1,1);
bar(x, p, 'k'); axis tight;
title('\rho_{0}')
subplot(2,1,2);
bar(y, q, 'k'); axis tight;
title('\rho_{1}')



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
cap = 1.3;
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
figure(2)
flux = plot(t,rk);
flux.LineWidth = 2;
flux.Color ='#0072BD';
title('Marginal of Time between [0,T]')

figure(3)
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
figure(4)
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

figure(5)
plot(t,rk)
%% Display the optimal assignement

figure(6)
p1 = plot3(xx1,xx2,p,'LineWidth',1.5,'Color','k');

hold on
p2 = plot3(yy1,yy2,q,'LineWidth',1.5,'Color','k');
tplot = plot3(tt2,tt1,rk,'LineWidth',1.5,'Color','k');
toll = plot3(tt2,tt1,tollheight,'LineWidth',1.5,'Color','r','LineStyle',':');
% bar3(tt2,tt1,tollheight);

% Plot the grey 'area' by Line
% pt Coupling (correct):
[Ipt,Jpt,~] = find(ptCoupling>0.0002);
for ipt=1:length(Ipt)
    pt_assignement = plot([xx(Ipt(ipt),1),tt2(Jpt(ipt))],[xx(Ipt(ipt),2),tt1(Jpt(ipt))],'LineWidth',3,'Color',[0.8 0.8 0.8]);
end

% qt Coupling (order fixed):
[Iqt,Jqt,~] = find(qtCoupling>0.0002);
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

% capacity = 2;
% yarea = [0.095 0.105 0.355 0.365 0.47 0.48 0.51 0.545 0.595 0.615 0.79 0.805];

% capacity = 1.2

% yarea = [0.03 0.045 0.865 0.875];
% xarea = zeros(size(yarea));
% zperiod = [0 capacity capacity 0];
% zarea = [zperiod];
% 
% tollgrey = fill3(xarea,yarea,zarea, [0.2 0.2 0.2]);
% tollgrey.FaceAlpha = 0.5;


ax = gca;
ax.FontSize = 13; 
hold off
xlabel('Space','FontSize',14)
ylabel('Time','FontSize',14)
% zlabel('Probability Density','FontSize',14)
zlabel('Density','FontSize',14)
view([-15 50])


%% Plot 3-D coupling/Sparse 3D matrix 
% Generate new coordinates for 3D plot
xx1 = x;
xx2 = ones(size(xx1));
xx = [xx1,xx2,p];
yy1 = y;
yy2 = zeros(size(yy1));
yy = [yy1,yy2,q];

tt1 = t;
tt2 = ones(size(tt1));
tollheight = capacity*ones(size(tt1));


%
list = find (Pi>0.00001);
[I,J,K] = ind2sub(size(Pi),list);
%
figure()
scalar = 9;
p1 = plot3(xx1,xx2,scalar.*p,'LineWidth',1.5,'Color','b');
hold on
p2 = plot3(yy2,yy1,scalar.*q,'LineWidth',1.5,'Color','r');
tplot = plot3(-scalar.*rk,tt2,tt1,'LineWidth',1.5,'Color','k');

 for i=1:length(I)
     % pt_assignement = plot3(x(I(i)),y(J(i)),t(K(i)),'.','Color','k','MarkerSize',4000*Pi(I(i),J(i),K(i)),'MarkerFaceColor','#D9FFFF');
     pt_assignement = plot3(x(I(i)),y(J(i)),t(K(i)),'.','MarkerSize',4000*Pi(I(i),J(i),K(i)),'MarkerFaceColor','#D9FFFF');
 end

hold off
xlabel('X','FontSize',14)
ylabel('Y','FontSize',14)
zlabel('T','FontSize',14)
view([-69 23])



















