%% Coupled Linear programming problem:
%%  1. Two cost matrices: C_{1}, C_{2}, needed to be defined properly
%%  2. Two coupling matrices: Pi_{1}, Pi_{2}
%%  3. Four marginals: p,q,u,v. Noticed that u,v is not given
%%  4. The upper bound r at z-axis
clear all
clc

%% Useful functions:
normalize = @(a)a/sum(a(:)); % Normalization function 
Gaussian = @(x,t0,sigma)exp(-(x-t0).^2/(2*sigma^2)); % Gaussian distribution generator

%% Define given the marginal disributions p,q in the problem
sigma0 = 0.01; sigma1 = .02; sigma2 = .03; sigma3 = .04; sigma4 = .05;sigma5= .08;

nx = 200;  % (not exceed 1000)
ny = 200;  % (not exceed 1000)
nz = 300;  % (not exceed 1000)

x = (1:nx)'/nx;  % x \in [0,1]
y = (1:ny)'/ny;  % y \in [0,1]
z = (1:nz)'/nz;  % z \in [0,1]
%z = (0:nz-1)'/nz+2; % z \in [2,3]
heightp = 1;
heightq = 5;

p = heightp.*Gaussian(y,.5,sigma3);
q = heightq.*Gaussian(x,.5,sigma5);

% Normalized!
p = normalize( p );
q = normalize( q );

figure(1)
subplot(2,1,1);
bar(x, p, 'k'); axis tight;
subplot(2,1,2);
bar(y, q, 'k'); axis tight;

%% Define the cost matrices C_{1} C_{2} in the problem 
[Z,X] = meshgrid(x,z);
[Z,Y] = meshgrid(y,z);
C1 = 1*(Z-X).^(2);
% C1 = 1*(Z.^2+X.^2);
% C1 = 1*(X.^2./Z);
% C1 = abs(Z.^(2)-X.^(2));

 C2 = 100*(Z-Y).^(2);
% C2 = 1*(Z.^2+Y.^2);
% C2 = 100*(Y.^2./Z);
% C2 = 20*abs(Z.^(2)-Y.^(2));
%% The upper bound on z-axis
r = 5;
r_unit = r/nz;

%% CVX solver
cvx_begin
cvx_solver mosek
    variable Pi_1(nx,nz)
    variable Pi_2(ny,nz)
    minimize( sum(sum(C1'.*Pi_1 + C2'.*Pi_2)) )
    subject to
     Pi_1*ones(nz,1) == p;
     Pi_2*ones(nz,1) == q;
     sum(Pi_1'*ones(nx,1)) == 1;
     sum(Pi_2'*ones(nx,1)) == 1;
     Pi_1 >= 0;
     Pi_2 >= 0;
     Pi_1'*ones(nx,1)+Pi_2'*ones(nx,1) <= r_unit;
cvx_end

%% Constraint check (suppose to be 0):
fprintf('Constraints deviation (should be 0):'); 
sum(Pi_1*ones(nz,1) - p)
% sum(Pi_1'*ones(nx,1)- u)
sum(Pi_2*ones(nz,1) - q)
% sum(Pi_2'*ones(ny,1)- v)


fprintf('Marginal at z-axis check (should be 1):'); 
sum(Pi_1'*ones(nx,1))
sum(Pi_2'*ones(nx,1))

%% Plotting the u,v
figure(2)
subplot(2,1,1);
bar(z, Pi_1'*ones(nx,1),'k'); % axis tight;
subplot(2,1,2);
bar(z, Pi_2'*ones(nx,1),'k'); % axis tight;
 
%% Plot in 3-D
scalar = 10;
xx = x;
xy = zeros(size(x));
xz = scalar.*p;

yx = zeros(size(y));
yy = y;
yz = scalar.*q;

zxp = scalar.* Pi_1'*ones(nx,1);
zxq = scalar.* Pi_2'*ones(nx,1);
zy = zeros(size(z));
zz = z;

figure(3)
p1 = plot3(xx,xy,xz,'LineWidth',2,'Color','r');
hold on
p2 = plot3(yx,yy,yz,'LineWidth',2,'Color','k');
p3 = plot3(zxp,zy,zz,'LineWidth',2,'Color','r');
p4 = plot3(zxq,zy,zz,'LineWidth',2,'Color','k');


% C1 = 2; C2=4
zareap = [0.3 0.3 0.4 0.4 0.6 0.6 0.7 0.7];
zareaq = [0.4 0.4 0.6 0.6];

yareap = zeros(size(zareap));
yareaq = zeros(size(zareaq));

xperiod = [0 scalar*r_unit scalar*r_unit 0];
xareap = [xperiod, xperiod];
xareaq = [xperiod];

areap = fill3(xareap,yareap,zareap,'k');
areap.FaceAlpha = 0.8;

areaq = fill3(xareaq,yareaq,zareaq, 'r');
areaq.FaceAlpha = 0.8;

ax = gca;
ax.FontSize = 13; 
hold off
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Coupled LP problem')
view([135 18])





















































