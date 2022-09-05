%% Multi-Marginal Optimal mass transport problem
% This script numerically solves an OMT problem as follows:
% 1. Transportation cost: C_{ijk} = \frac{x_{i}^2}{t_{k}} + \frac{y_{i}^2}{T - t_{k}}
% 2. Transportation map: Pi_{ijk} with its three marginials: p_i, q_i, r_k
%    such that \sum_{jk} = p_{i}, \sum_{ik} = q_{j}, \sum_{ij} = r_k.
%    p_i and q_j are two probability vector   
% 3. The mariginal along the t cooridinates has an upper-bound: \frac{1.1}{T}.
% 4. p is transported to q in the time interval [0,T]


% The optimizer is the vector r_k which minimize the transportation cost:
% min_{ijk} C_{ijk}.*\Pi_{ijk}

% CVX with MOSEK solver is used:
% Installation of CVX: http://cvxr.com/cvx/

clear all
clc

%% Define the two probability distributions p_{i} and q_{j} and cooridinates X(p), Y(q)
% Consider p_{i} and q_{j} are two discretes distributions 
% p_{i}: i = 1, 2,..., n1
% q_{j}: j = 1, 2,..., n2

%% Gaussian distributions p,q in histogram

% 
n1 = 200;  % (not exceed 1000)
n2 = 200;  % (not exceed 1000)

x = (0:n1-1)'/n1-1; % x \in [-1,0]
y = (1:n2)'/n2;     % y \in [ 0,1]

% Useful functions 
normalize = @(a)a/sum(a(:));                           % Normalization function 
Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) ); % Gaussian distribution generator

sigma1 = .03; sigma2 = .04; sigma3 = .01;
% Marginal p
% p = Gaussian(x,-.9,sigma1);
p = Gaussian(x, -.8,sigma1) + Gaussian(x, -.5,sigma2);
% Marginal q
q = Gaussian(y, .4,sigma2);
% q = Gaussian(y, .2,sigma2) + Gaussian(y, .8,sigma3);

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
% The capacity constraint should be larger than 1

% Define the capacity 
cap = 0.55;

cvx_begin
cvx_solver mosek
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
      Pi >= 0;
      0 <= rk <= cap*T/nt;
cvx_end


%% Check 3 marginal constraints:
fprintf('Values should be 0:'); 
sum(squeeze(sum(sum(permute(Pi,[2,3,1])))) - p)
sum(squeeze(sum(sum(permute(Pi,[3,1,2])))) - q)
sum(squeeze(sum(sum(Pi)))- rk)

%% Show that t marginal is also a measure
fprintf('Values should be 1:'); 
sum(squeeze(sum(sum(Pi))))

%% Plots of rk, evolutionary probability distributions at moment t
figure(4)
flux = plot(t,rk);
flux.LineWidth = 2;
flux.Color ='#0072BD';
title('Marginal of Time between [0,T]')


figure(5)
bar(t,rk, 'k'); axis tight;
title('Marginal of Time between [0,T]')


%% 3 Marginals Plot




%% Wasserstein-2 distance Computation and Comparision





%% velocity field plot 
























