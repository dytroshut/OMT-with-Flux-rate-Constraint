%% Generalized multimarginal problem:
%% Consider the Generalized multimarginal problem, in which we now have two 
%% (or even more) toll placed at location, say x_{0} and y_{0}. The time 
%% variable now becomes t_{1} and t_{2} which is the moments of the particles 
%% throught the first and the second toll. The cost matrix is defined as follows:
% C_{x,y,t1,t2} = (x-x_{0})^{2}/t_1 + (y_{0}-x_{0})^{2}/t_2 +
% (y-y_{0})^{2}/(T-t_1-t_2)

% In the experiment, we set T=1, x_{0} = -0.2, y_{0} = 0.2;

clear all
clc


n1 = 40;  % (not exceed 1000)
n2 = 40;  % (not exceed 1000)

x = (0:n1-1)'/n1-1; % x \in [-1,0]
y = (1:n2)'/n2;     % y \in [ 0,1]

normalize = @(a)a/sum(a(:));  % Normalization function 

% Generate Dirac
p = zeros(size(x));
p((n1/5-2):(n1/5+2)) = 2;      % 'Dirac 1'

q = zeros(size(y));

% q((n2*3/5-2):(n2*3/5+2)) = 2;      % 'Dirac 2'
q((n2*4/5-2):(n2*4/5+2)) = 2;  % 'Dirac 3'

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
nt1 = 40;
nt2 = 40;
% Define the time interval [0,T], notice that T is not neccesary to be 1.
T = 1;
t1 = 1/2*T*(0+0.01:nt1-0.01)'/nt1;
t2 = 1/2*T*(0+0.01:nt2-0.01)'/nt2+0.5;

%% Define the cost function C(x,y,t):
% C(x,y,t) = C_{ijk} = (x_i)^2./t_k + (y_i)^2/(T - t_k)
% Compute the 3-dimensional Cost matrix
[X,Y,T1,T2] = ndgrid(x,y,t1,t2);
C = (X+0.4).^2./T1 + (0.8).^2./(T2-T1) + (Y-0.4).^2./(T-T2);
% C = (X+0.2).^2./T1 + (0.4).^2./(T2-T1)  +(Y-0.2).^2./(T-T2);
% C = (X+0.2).^2./T1 + (0.4).^2.*inv_pos(T2-T1) +(Y-0.2).^2./(T-T2);
%% CVX solver 
% The CVX solver: MOSEK needs to be installed. 
% The sum function in the CVX package do not support high dimension summation
% The capacity constraint should be larger than 1: cap*T/nt>1

% Define the capacity 
cap = 1.5;
capacity1 = cap*T/nt1;
capacity2 = cap*T/nt2;

cvx_begin
cvx_solver mosek
cvx_precision best  % Options: low, medium, default, high, best
    variable Pi(n1,n2,nt1,nt2)
    variable rk1(nt1)
    variable rk2(nt2)
    minimize( sum(sum(sum(sum(C.*Pi)))) )
    subject to  
    % Define multi-marginal coupling counstraints:
    % 1.Marginal constraints on x coordinate
      squeeze( sum(sum(sum(permute(Pi,[2,3,4,1])))) ) == p;
    % 2.Marginal constraints on y coordinate   
      squeeze( sum(sum(sum(permute(Pi,[3,4,1,2])))) ) == q;
    % 3.Marinigal on t1-coordinate
      squeeze( sum(sum(sum(permute(Pi,[4,1,2,3])))) ) == rk1;
    % 4.Marinigal on t1-coordinate
      squeeze( sum(sum(sum(Pi))) ) == rk2; 
    % 4.Element-wise lower bound for \Pi
      Pi >= 0;
    % 5.Capacity/Toll constraints
        0 <= rk1 <= capacity1;    % r1_{k} with upperbound,
        %0 <= rk1;                  % r1_{k} without upperbound.
        0 <= rk2 <= capacity2;    % r2_{k} with upperbound,
        %0 <= rk2;                  % r2_{k} without upperbound.      
cvx_end


%% Constraints and Marginals checking
fprintf('Constraints deviation (should be 0):'); 
sum( squeeze( sum(sum(sum(permute(Pi,[2,3,4,1])))) ) - p )
sum( squeeze( sum(sum(sum(permute(Pi,[3,4,1,2])))) ) - q )
sum( squeeze( sum(sum(sum(permute(Pi,[4,1,2,3])))) )- rk1 )
sum( squeeze( sum(sum(sum(Pi))) ) - rk2 )

% Show that t marginal is a prbability measure
fprintf('Time marginals probability vector check (should be 1):'); 
sum(rk1)
sum(rk2)

%% testing t-marginals

figure()
subplot(2,1,1);
bar(t1,rk1,'k');
title('r1_{t}')
subplot(2,1,2);
bar(t2,rk2,'k');
title('r2_{t}')


%% Generate new coordinates for 3D plot
xx = x;
xy = zeros(size(xx));
xp = [xx,xy,p];

yy = y;
yx= ones(size(yy));
yq = [yy,yx,q];

toll_1_loc = -0.4;
toll_2_loc = 0.4;
ty1 = t1;
tx1 = (toll_1_loc).*ones(size(ty1));

ty2 = t2;
tx2 = (toll_2_loc).*ones(size(ty2));

%% figure: Line plot3
figure()
p1 = plot3(xx,xy,p,'LineWidth',2);
hold on
p2 = plot3(yy,yx,q,'LineWidth',2);
toll_1 = plot3(tx1,ty1,rk1,'LineWidth',1,'Color','k');
toll_2 = plot3(tx2,ty2,rk2,'LineWidth',1,'Color','k');

hold off
grid on
xlabel('Space')
ylabel('Time')
zlabel('Probability Density')
title('Two tolls')

%% Couplings 
pt1Coupling = squeeze(sum(Pi,[2,4]));
qt2Coupling = squeeze(sum(Pi,[1,3]));
tollCoupling = squeeze(sum(Pi,[1,2]));

figure()
p1 = plot3(xx,xy,p,'LineWidth',2,'Color','k');
hold on
p2 = plot3(yy,yx,q,'LineWidth',2,'Color','k');

toll_1 = plot3(tx1,ty1,rk1,'LineWidth',1,'Color','k');
toll_2 = plot3(tx2,ty2,rk2,'LineWidth',1,'Color','k');


[Ipt,Jp,~] = find(pt1Coupling);
for ipt=1:length(Ipt)
    pt_assignement = plot([xp(Ipt(ipt),1),tx1(Jp(ipt))],[xp(Ipt(ipt),2),ty1(Jp(ipt))],'LineWidth',3,'Color',[0.8 0.8 0.8]);
end

% qt Coupling (order fixed):
[Iqt,Jqt,~] = find(qt2Coupling);
for iqt=1:length(Iqt)
    qt_assignement = plot([yq(Iqt(iqt),1),tx2(Jqt(iqt))],[yq(Iqt(iqt),2),ty2(Jqt(iqt))],'LineWidth',3,'Color',[0.8 0.8 0.8]);
end

% tolls Coupling (order fixed):
[Itt,Jtt,~] = find(tollCoupling);
for iqt=1:length(Itt)
    qt_assignement = plot([tx1(Itt(iqt)),tx2(Jtt(iqt))],[ty1(Itt(iqt)),ty2(Jtt(iqt))],'LineWidth',3,'Color',[0.8 0.8 0.8]);
end

% t-axis plot
plot([toll_1_loc,toll_1_loc],[0,1],'LineWidth',2,'Color','k');
plot([toll_2_loc,toll_2_loc],[0,1],'LineWidth',2,'Color','k');

yarea1 = [.11 .11 .43 .43];
yarea2 = [.56 .56 .8876 .8876];
xarea1 = toll_1_loc*ones(size(yarea1));
xarea2 = toll_2_loc*ones(size(yarea2));
zarea = [0 capacity1 capacity1 0];

toll_1 = fill3(xarea1,yarea1,zarea, [0.6350 0.0780 0.1840]);
toll_1.FaceAlpha = 0.7;
toll_2 = fill3(xarea2,yarea2,zarea, [0 0.4470 0.7410]);
toll_2.FaceAlpha = 0.7;

ax = gca;
ax.FontSize = 13; 
hold off
xlabel('Space','FontSize',14)
ylabel('Time','FontSize',14)
zlabel('Density','FontSize',14)
view([-15 50])











