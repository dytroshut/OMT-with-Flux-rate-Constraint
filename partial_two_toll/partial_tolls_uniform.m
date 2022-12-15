%% Two toll with one separating masses: every time interval is [0,1]

%% Two tolls with one separating masses

%% We consider the OMT problem with two tolls but with one separating mass, 
% which says that....


%% Define p,q,u,v and 
clear all
clc

% we use the uniform distribution for the sake of simplification:
n1 = 50;  % (not exceed 1000)
n2 = 50;  % (not exceed 1000)

x = (0:n1-1)'/n1-1; % x \in [-1,0]
y = (1:n2)'/n2;     % y \in [ 0,1]

% Useful functions 
normalize = @(a)(0.5)*a/sum(a(:));  % Normalization function 


% Generate Dirac
p = zeros(size(x));
p((n1/5-3):(n1/5+3)) = 2;      % 'p'

q = zeros(size(x));
q((n1*4/5-3):(n1*4/5+3)) = 2;  % 'q'

v = zeros(size(y));
v((n2/2+10):(n2/2+18)) = 2;      % 'Dirac 2'

u = zeros(size(y));
u((n2/2-14):(n2/2-6)) = 2;      % 'Dirac 2'

p = normalize(p);
q = normalize(q);
u = normalize(u);
v = normalize(v);

%%
% Plot the two distributions p and q
figure(1)
subplot(2,1,1);
hold on
bar(x, p,'k'); 
bar(x, q,'r');
hold off
axis tight;
title('p, q')
subplot(2,1,2);
hold on
bar(y, u, 'k'); 
bar(y, v, 'r');
hold off
axis tight;
title('u, v')


%% Define the time t_{i} \in [0, T] (tk):
% Discretization of the time mariginal by nt (not exceed 500)
nt1 = 50;
nt2 = 50;
nt = 50; 
% Be careful with nt1, nt2, nt, dimension need to be matched nt = 2*nt1 =2*nt2 

% Define the time interval [0,T], notice that T is not neccesary to be 1.
T = 1;
t1 = T*(1-0.01:nt1-0.01)'/nt1;
t2 = T*(1-0.015:nt2-0.015)'/nt2;
t = T*(0+0.01:nt-0.01)'/nt;

%% Define the weighted cost function C(x,y,t):
% C(x,y,t) = C_{ijk} = (x_i)^2./t_k + (y_i)^2/(T - t_k)
% Compute the 3-dimensional Cost matrix
[X,Y,T1,T2] = ndgrid(x,y,t1,t2);

Cpq = (X+0.5).^2./T1 + (Y).^2./(T-T2);
Ctt = (0.5).^2./(T2-T1);

Ctt(Ctt<=0) = 10e4;
C1= Cpq+Ctt;

C1(C1==Inf) = 10e4;
C1(isnan(C1)) = 0;

[X,Y,tk] = meshgrid(x,y,t);
C2 = X.^2./tk + Y.^2./(T-tk);


%% CVX solver 
% The CVX solver: MOSEK needs to be installed. 
% The sum function in the CVX package do not support high dimension summation
% The capacity constraint should be larger than 1: cap*T/nt>1

% Define the capacity 
cap1 = 1.2;
cap2 = 1.2;
capacity1 = cap1*T/nt1;
capacity2 = cap2*T/nt;

cvx_begin
cvx_solver mosek
cvx_precision best  % Options: low, medium, default, high, best
    variable Pi1(n1,n2,nt1,nt2)
    variable Pi2(n1,n2,nt)
    variable rk1(nt1)
    variable rk2(nt)
    minimize( sum( squeeze(sum(sum(sum(C1.*Pi1)))) )+ sum( squeeze(sum(sum(C2.*Pi2))) ) )
    subject to  
    % 1.Marginal constraints on x coordinate Pi1
      squeeze( sum(sum(sum(permute(Pi1,[2,3,4,1])))) ) == p;
    % 2.Marginal constraints on y coordinate Pi1 
      squeeze( sum(sum(sum(permute(Pi1,[3,4,1,2])))) ) == u;
    % 3.Marginal constraints on x coordinate Pi2
      squeeze(sum(sum(permute(Pi2,[2,3,1])))) == q;
    % 4.Marginal constraints on y coordinate Pi2       
      squeeze(sum(sum(permute(Pi2,[3,1,2])))) == v;         
    % 5.Marinigal on t1-coordinate
      squeeze( sum(sum(sum(permute(Pi1,[4,1,2,3])))) ) == rk1; 
      
    % 6.Marinigal on t2-coordinate
      squeeze( sum(sum(sum(Pi1)))) + squeeze(sum(sum(Pi2))) == rk2; 
    % 7.Element-wise lower bound for \Pi
      Pi1 >= 0;
      Pi2 >= 0;
    % 8.Capacity/Toll constraints
       %0 <= rk1 <= capacity1;    % r1_{k} with upperbound,
       %0 <= rk2 <= capacity2;    % r2_{k} with upperbound,
       0 <= rk1;                  % r1_{k} without upperbound.
       0 <= rk2;                  % r2_{k} without upperbound.      
cvx_end



%% Constraints and Marginals checking
fprintf('Constraints deviation (should be 0):'); 
sum( squeeze( sum(sum(sum(permute(Pi1,[2,3,4,1])))) ) - p )
sum( squeeze( sum(sum(sum(permute(Pi1,[3,4,1,2])))) ) - u )
sum( squeeze( sum(sum(sum(permute(Pi1,[4,1,2,3])))) )- rk1 )

sum( squeeze(sum(sum(permute(Pi2,[2,3,1])))) - p )
sum( squeeze(sum(sum(permute(Pi2,[3,1,2])))) - u )

% Show that t marginal is a prbability measure
fprintf('Time marginals probability vector check (should be 1/2):'); 
sum(rk1)
fprintf('Time marginals probability vector check (should be 1 ):'); 
sum(rk2)


%% testing t-marginals

figure()
subplot(2,1,1);
bar(t1,rk1,'k');
title('r1_{t}')
subplot(2,1,2);
bar(t,rk2,'k');
title('r2_{t}')


%% Generate new coordinates for 3D plot
xx = x;
xy = zeros(size(xx));
xp = [xx,xy,p];
xq = [xx,xy,q];
yy = y;
yx= ones(size(yy));
yu = [yy,yx,u];
yv = [yy,yx,v];
toll_1_loc = -0.5;
toll_2_loc = 0;
ty1 = t1;
tx1 = (toll_1_loc).*ones(size(ty1));
ty2 = t2;
tx2 = (toll_2_loc).*ones(size(ty2));
ty =  t;
tx = (toll_2_loc).*ones(size(ty));
rk21 = squeeze( sum(sum(sum(Pi1))) );
rk22 = squeeze(sum(sum(Pi2)));

%% Grey and Pink trajectory plot
pt1Coupling = squeeze(sum(Pi1,[2,4]));
ut2Coupling = squeeze(sum(Pi1,[1,3]));
tollCoupling = squeeze(sum(Pi1,[1,2]));

qtCoupling = squeeze(sum(Pi2,2));
vtCoupling = squeeze(sum(Pi2,1));

figure()
p1 = plot3(xx,xy,p,'LineWidth',2,'Color','k');
hold on
p2 = plot3(yy,yx,u,'LineWidth',2,'Color','k');
p3 = plot3(xx,xy,q,'LineWidth',2,'Color','#A2142F');
p4 = plot3(yy,yx,v,'LineWidth',2,'Color','#A2142F');

toll_1 = plot3(tx1,ty1,rk1,'LineWidth',2,'Color','k');
toll_21 = plot3(tx2,ty2,rk21,'LineWidth',2,'Color','k');
toll_22 = plot3(tx,ty,rk22,'LineWidth',2,'Color','#A2142F');


% pt Coupling (order fixed):
[Ipt,Jp,~] = find(pt1Coupling);
for ipt=1:length(Ipt)
    pt_assignement = plot([xp(Ipt(ipt),1),tx1(Jp(ipt))],[xp(Ipt(ipt),2),ty1(Jp(ipt))],'LineWidth',3,'Color',[0.6 0.6 0.6]);
end

% ut Coupling (order fixed):
[Iut,Jut,~] = find(ut2Coupling);
for iut=1:length(Iut)
    ut_assignement = plot([yu(Iut(iut),1),tx2(Jut(iut))],[yu(Iut(iut),2),ty2(Jut(iut))],'LineWidth',3,'Color',[0.6 0.6 0.6]);
end

% tolls Coupling (order fixed):
[Itt,Jtt,~] = find(tollCoupling);
for itt=1:length(Itt)
    tt_assignement = plot([tx1(Itt(itt)),tx2(Jtt(itt))],[ty1(Itt(itt)),ty2(Jtt(itt))],'LineWidth',3,'Color',[0.6 0.6 0.6]);
end

% qt Coupling
[Iqt,Jqt,~] = find(qtCoupling);
for iqt=1:length(Iqt)
    qt_assignement = plot([xq(Iqt(iqt),1),tx(Jqt(iqt))],[xq(Iqt(iqt),2),ty(Jqt(iqt))],'LineWidth',3,'Color','#FFCCCB');
end

% vt Coupling (order fixed):
[Ivt,Jvt,~] = find(vtCoupling);
for ivt=1:length(Ivt)
    vt_assignement = plot([yv(Ivt(ivt),1),tx(Jvt(ivt))],[yv(Ivt(ivt),2),ty(Jvt(ivt))],'LineWidth',3,'Color','#FFCCCB');
end



hold off
xlabel('Space','FontSize',14)
ylabel('Time','FontSize',14)
zlabel('Density','FontSize',14)
view([-4 41])