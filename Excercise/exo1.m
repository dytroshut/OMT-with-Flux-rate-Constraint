clear all
clc

%% Useful functions:
normalize = @(a)a/sum(a(:)); % Normalization function 
Gaussian = @(x,t0,sigma)exp(-(x-t0).^2/(2*sigma^2)); % Gaussian distribution generator

%% Define given the marginal disributions p,q in the problem
sigma0 = 0.01; sigma1 = .02; sigma2 = .03; sigma3 = .04; sigma4 = .05;sigma5= .06;

nx = 200;  % (not exceed 1000)
ny = 200;  % (not exceed 1000)

x = (1:nx)'/nx;  % x \in [0,1]
y = (1:ny)'/ny;  % y \in [0,1]

heightp = 1;
p = heightp.*Gaussian(y,.5,0.2);

% Normalized!
p = normalize(p);


figure(1)
bar(x, p, 'k'); axis tight;


%% Define the cost matrices C_{1} C_{2} in the problem 
[Y,X] = meshgrid(x,y);
% C = abs((Y-X).^(3));
C = (Y-X).^(2);
% C = X.^2./Y;
%% The upper bound on z-axis
r = 1.6;
r_unit = r/ny;

%% CVX solver
cvx_begin
    variable Pi(nx,ny)
    minimize( sum(sum(C'.*Pi)) )
    subject to
     Pi*ones(ny,1) == p;
     sum(Pi'*ones(nx,1)) == 1;
     Pi'*ones(nx,1) <= r_unit;
     Pi>=0;
cvx_end

%% Constraint check (suppose to be 0):
fprintf('Constraints deviation (should be 0):'); 
sum(Pi*ones(ny,1) - p)


fprintf('Marginal at z-axis check (should be 1):'); 
sum(Pi'*ones(nx,1))

%% Plotting the u,v
figure()
plot(y, Pi'*ones(nx,1)); axis tight;


 
