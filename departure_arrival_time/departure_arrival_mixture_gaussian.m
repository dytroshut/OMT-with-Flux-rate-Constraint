%% Relaxed Gaussian transportation
%% OMT: Departure and Arriving at different time
% In this generalization?we consider the case when the particles departure 
% and arrival at different time. Different from the classic optimal mass
% transport problem, in which all the particles have to leave at t0 and
% arrive at tf, we consider a constraint so that the their is limitation of
% how much the mass can leave/arrive in a unit of time. The code provides
% a numerical approach of such problem so that we consider (x,y,t0,t1) so 
% that x,y is the location of the particles and t0(x) is the departure time
% of mass at x and t1(y) is the arrival time of the particle goes to y.

% Moreover, t0,t1 can be bounded so that instead of giving the explict
% departure/arrival times, we put a constraint of the quatity of the particle
% that allowed to leave/arrive at the time interval.

%% Define the p,q at x,y and the time function t0,t1
clear all
clc

n1 = 30;  % (not exceed 1000)
n2 = 30;  % (not exceed 1000)

x = (0.0001:n1-1+0.0001)'/n1-1; % x \in [-1,0]
y = (1-0.0002:n2-0.0002)'/n2;     % y \in [ 0,1]

% Useful functions 
normalize = @(a)a/sum(a(:));                           % Normalization function 
Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) ); % Gaussian distribution generator

sigma1 = .02; sigma2 = .03; sigma3 = .04; sigma4 = .1;
% Marginal p
p = 2.*Gaussian(x, -.4, 0.08) + 3.*Gaussian(x, -.6,.05 );
% Marginal q
q = 3.*Gaussian(y, .4,  0.1) + 5.*Gaussian(y, .7, .05);



p = normalize(p);
q = normalize(q);

% Plot the two distributions p and q
figure()
subplot(2,1,1);
bar(x, p, 'k'); axis tight;
title('\rho_{0}')
subplot(2,1,2);
bar(y, q, 'k'); axis tight;
title('\rho_{1}')


%% Define the time t_{i} \in [0, T] (tk):
% Discretization of the time mariginal by nt (not exceed 500)
nt0 = 100;
nt1 = 100;
% Define the time interval [0,T], notice that T is not neccesary to be 1.
T = 1;
% t0 = 1/2*T*(0+0.01:nt0-0.01)'/nt0;
% t1 = 1/2*T*(0+0.01:nt1-0.01)'/nt1+0.5;

t0 = T*(1-0.01:nt0-0.01)'/nt0;
t1 = T*(1-0.01:nt1-0.01)'/nt1;
%% Define the weighted cost function C(x,y,t):

%% C(x,y,t0,t1) = (y-x)^2/(t1-t0) + extra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% epsilon = 1;
[X,Y,T0,T1] = ndgrid(x,y,t0,t1);

C = T0./X.^2 + (Y-X).^2 - T1./Y.^2;
% C = -exp(X).*exp(-T0) + (Y-X).^2 + exp(-Y).*exp(-T1);

%% CVX solver 
% The CVX solver: MOSEK needs to be installed. 
% The sum function in the CVX package do not support high dimension summation
% The capacity constraint should be larger than 1: cap*T/nt>1

% Define the capacity 
cap1 = 1.1;  % nice bound cap = 1.6 ,3
cap2 = 5;  % nice bound cap = 1.6 ,3
capacity1 = cap1*T/nt0;
capacity2 = cap2*T/nt1;

cvx_begin
cvx_solver mosek
cvx_precision best  % Options: low, medium, default, high, best
    variable Pi(n1,n2,nt0,nt1)
    variable Dt(nt0)
    variable At(nt1)
    minimize( sum( sum( sum( sum(C.*Pi) ) ) ) )
    subject to  
    % Define multi-marginal coupling counstraints:
    % 1.Marginal constraints on x coordinate
      squeeze( sum(sum(sum(permute(Pi,[2,3,4,1])))) ) == p;
    % 2.Marginal constraints on y coordinate   
      squeeze( sum(sum(sum(permute(Pi,[3,4,1,2])))) ) == q;
    % 3.Marinigal on t0-coordinate
      squeeze( sum(sum(sum(permute(Pi,[4,1,2,3])))) ) == Dt;
    % 4.Marinigal on t1-coordinate
      squeeze( sum(sum(sum(Pi))) ) == At; 
    % 4.Element-wise lower bound for \Pi
      Pi >= 0;
    % 5.Capacity/Toll constraints
      0 <= Dt <= capacity1;    % r1_{k} with upperbound,
      % 0 <= Dt;                  % r1_{k} without upperbound.
      0 <= At <= capacity2;    % r2_{k} with upperbound,
      % 0 <= At;                  % r2_{k} without upperbound.      
cvx_end


%% Constraints and Marginals checking
fprintf('Constraints deviation (should be 0):'); 
sum( squeeze( sum(sum(sum(permute(Pi,[2,3,4,1])))) ) - p )
sum( squeeze( sum(sum(sum(permute(Pi,[3,4,1,2])))) ) - q )
sum( squeeze( sum(sum(sum(permute(Pi,[4,1,2,3])))) )- Dt )
sum( squeeze( sum(sum(sum(Pi))) ) - At )

% Show that t marginal is a prbability measure
fprintf('Time marginals probability vector check (should be 1):'); 
sum(Dt)
sum(At)


%% Plot x-t0 coupling and y-t1 coupling
pt0Coupling = squeeze(sum(Pi,[2,4]));
qt1Coupling = squeeze(sum(Pi,[1,3]));
pqCoupling  = squeeze(sum(Pi,[3,4]));


figure()
hold on 
area(x, -p, 'FaceColor', [0 0.4470 0.7410]);
plot(Dt, t0,'LineWidth',2, 'Color', [0.4660 0.6740 0.1880]);
% area(Dt',t0');
XX1 = [0 max(Dt) max(Dt) 0];
list1 = find(Dt==max(Dt));
list2 = find(Dt==0);
YY1 = [t0(list1(1)) t0(list1(1)) t0(list1(end)) t0(list1(end))];
patch(XX1,YY1,[0.4660 0.6740 0.1880]);

% Coupling plot
[I,J] = find(pt0Coupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(x(I(i)),t0(J(i)),'o','Color','k','MarkerSize',...
     200*pt0Coupling(I(i),J(i)),'MarkerFaceColor','k');
 end

hold off
xlabel('X','FontSize',18)
ylabel('T_{dep}','FontSize',18)
axis tight
title('Coupling: (\mu, t_{dep})','FontSize',18)

%%
figure()
hold on
area(y, -q, 'FaceColor', [0.6350 0.0780 0.1840]);
plot(-At,t1,'LineWidth',2,'Color', [0.4660 0.6740 0.1880]);

XX2 = 1*[0 -max(At) -max(At) 0];
list1 = find(At==max(At));
list2 = find(At==0);
YY2 = [t1(list1(1)) t1(list1(1)) t1(list1(end)) t1(list1(end))];
patch(XX2,YY2,[0.4660 0.6740 0.1880]);

[I,J] = find(qt1Coupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(y(I(i)),t1(J(i)),'o','Color','k','MarkerSize',...
     300*qt1Coupling(I(i),J(i)),'MarkerFaceColor','k');
 end
hold off
axis tight;
title('Coupling: (\nu, t_{arr})','FontSize',18)
xlabel('Y','FontSize',18)
ylabel('T_1','FontSize',18)


%%
figure()
hold on
plot(x, p, 'LineWidth',2);
plot(-q, y, 'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);

[I,J] = find(pqCoupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(x(I(i)),y(J(i)),'-o','Color','k','MarkerSize',...
     200*pqCoupling(I(i),J(i)),'MarkerFaceColor','k');
 end
 
hold off
axis tight;
title('Coupling:(\mu,\nu)','FontSize',18)
xlabel('X','FontSize',18)
ylabel('Y','FontSize',18)

% view(-90, 90) %# Swap the axes
% set(gca, 'ydir', 'reverse'); %# Reverse the y-axis (Optional step)



%% One figure for paper
figure()
subplot(3,3,[1,5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot(x, p, 'LineWidth',2);
plot(-q, y, 'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);

[I,J] = find(pqCoupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(x(I(i)),y(J(i)),'-o','Color','k','MarkerSize',...
     100*pqCoupling(I(i),J(i)),'MarkerFaceColor','k');
 end
 
hold off
axis tight;
ttl = title('(\mu,\nu)','FontSize',12);
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)



subplot(3,3,[7,8]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on 
area(x, -p, 'FaceColor', [0 0.4470 0.7410]);
plot(Dt, t0,'LineWidth',2, 'Color', [0.4660 0.6740 0.1880]);
% area(Dt',t0');
XX1 = 3*[0 max(Dt) max(Dt) 0];
list1 = find(Dt==max(Dt));
list2 = find(Dt==0);
YY1 = [t0(list1(1)) t0(list1(1)) t0(list1(end)) t0(list1(end))];
patch(XX1,YY1,[0.4660 0.6740 0.1880]);

% Coupling plot
[I,J] = find(pt0Coupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(x(I(i)),t0(J(i)),'o','Color','k','MarkerSize',...
     120*pt0Coupling(I(i),J(i)),'MarkerFaceColor','k');
 end

hold off
xlabel('x','FontSize',12)
ylabel('t_{dep}','FontSize',12)
axis tight
%title('(\mu, t_{dep})','FontSize',12)



subplot(3,3,[3,6]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
area(y, -q, 'FaceColor', [0.6350 0.0780 0.1840]);
plot(-At,t1,'LineWidth',2,'Color', [0.4660 0.6740 0.1880]);

XX2 = 3*[0 -max(At) -max(At) 0];
list1 = find(At==max(At));
list2 = find(At==0);
YY2 = [t1(list1(1)) t1(list1(1)) t1(list1(end)) t1(list1(end))];
patch(XX2,YY2,[0.4660 0.6740 0.1880]);

[I,J] = find(qt1Coupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(y(I(i)),t1(J(i)),'o','Color','k','MarkerSize',...
     80*qt1Coupling(I(i),J(i)),'MarkerFaceColor','k'); % 200 radius for plot 1.6
 end
 
hold off
view(90, -90) %# Swap the axes
axis tight;
%title('(\nu, t_{arr})','FontSize',12)
xlabel('y','FontSize',12)
ylabel('t_{arr}','FontSize',12)






%% One figure for paper (version 2)
figure()
subplot(3,3,[1,5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot(x, p, 'LineWidth',2);
plot(-q, y, 'LineWidth',2,'Color',[0.6350 0.0780 0.1840]);

[I,J] = find(pqCoupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(x(I(i)),y(J(i)),'-o','Color','k','MarkerSize',...
     100*pqCoupling(I(i),J(i)),'MarkerFaceColor','k');
 end
 
hold off
axis tight;
% xlabel('x','FontSize',12)
ylabel('y','FontSize',12)



subplot(3,3,[7,8]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on 
area(x, p, 'FaceColor', [0 0.4470 0.7410]);
plot(-Dt, t0,'LineWidth',2, 'Color', [0.4660 0.6740 0.1880]);
% area(Dt',t0');
XX1 = -3*[0 max(Dt) max(Dt) 0];
list1 = find(Dt==max(Dt));
list2 = find(Dt==0);
YY1 = [t0(list1(1)) t0(list1(1)) t0(list1(end)) t0(list1(end))];
patch(XX1,YY1,[0.4660 0.6740 0.1880]);

% Coupling plot
[I,J] = find(pt0Coupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(x(I(i)),t0(J(i)),'o','Color','k','MarkerSize',...
     120*pt0Coupling(I(i),J(i)),'MarkerFaceColor','k');
 end

hold off
xlabel('x','FontSize',12)
ylabel('t_{d}','FontSize',12)
axis tight


subplot(3,3,[3,6]);
hold on
area(y, q, 'FaceColor', [0.6350 0.0780 0.1840]);
plot(At,t1,'LineWidth',2,'Color', [0.4660 0.6740 0.1880]);

XX2 = -3*[0 -max(At) -max(At) 0];
list1 = find(At==max(At));
list2 = find(At==0);
YY2 = [t1(list1(1)) t1(list1(1)) t1(list1(end)) t1(list1(end))];
patch(XX2,YY2,[0.4660 0.6740 0.1880]);

[I,J] = find(qt1Coupling>0.005);
 for i=1:length(I)
     pt_assignement = plot(y(I(i)),t1(J(i)),'o','Color','k','MarkerSize',...
     80*qt1Coupling(I(i),J(i)),'MarkerFaceColor','k'); % 200 radius for plot 1.6
 end
 
hold off
view(90, -90)
axis tight;
%xlabel('y','FontSize',12)
ylabel('t_{a}','FontSize',12)




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure()
% hold on 
% area(x,-p);
% plot(Dt,t0,'LineWidth',2);
% % area(Dt',t0');
% XX1 = [0 max(Dt) max(Dt) 0];
% list1 = find(Dt==max(Dt));
% list2 = find(Dt==0);
% % YY1 = [0 t0(list1(1)) t0(list1(end)) t0(list2(1))];
% % patch(XX1,YY1,[0.8500 0.3250 0.0980]);
% 
% % Coupling plot
% [I,J] = find(pt0Coupling>0.005);
%  for i=1:length(I)
%      pt_assignement = plot(x(I(i)),t0(J(i)),'o','Color','k','MarkerFaceColor','k');
%  end
% 
% hold off
% xlabel('X','FontSize',18)
% ylabel('T_{dep}','FontSize',18)
% axis tight
% title('Coupling: (\mu, t_{dep})','FontSize',18)

%%

% figure()
% hold on
% area(y, -q);
% plot(-At,t1,'LineWidth',2);
% 
% XX2 = [0 -max(At) -max(At) 0];
% list1 = find(At==max(At));
% list2 = find(At==0);
% % YY2 = [t1(list2(end-1)) t1(list1(1)) t1(list1(end)) 1];
% % patch(XX2,YY2,[0.8500 0.3250 0.0980]);
% 
% [I,J] = find(qt1Coupling>0.005);
%  for i=1:length(I)
%      pt_assignement = plot(y(I(i)),t1(J(i)),'o','Color','k','MarkerFaceColor','k');
%  end
% hold off
% axis tight;
% title('Coupling: (\nu, t_{arr})','FontSize',18)
% xlabel('Y','FontSize',18)
% ylabel('T_1','FontSize',18)


%%

% figure()
% hold on
% plot(x, p, 'LineWidth',2);
% plot(-q, y, 'LineWidth',2);
% [I,J] = find(pqCoupling>0.005);
%  for i=1:length(I)
%      pt_assignement = plot(x(I(i)),y(J(i)),'-o','Color','k','MarkerFaceColor','k');
%  end
%  
% hold off
% axis tight;
% title('Coupling:(\mu,\nu)','FontSize',18)
% xlabel('X','FontSize',18)
% ylabel('Y','FontSize',18)
















