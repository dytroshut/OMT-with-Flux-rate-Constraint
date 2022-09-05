%% Linear programming problem/ Coupled LP problem

% This script numerically solves an OMT problem as follows:
% 1. Transportation cost: C_{ij} 
% 2. Transportation map: Pi_{ij} with its two marginials: p_i, q_i
%    such that \sum_{j} = p_{i}, \sum_{i} = q_{j}
%    p_i and q_j are two probability vector   
% 3. The mariginal along the t cooridinates has an upper-bound

% The optimizer is the vector r_k which minimize the transportation cost:
% min_{ij} C_{ij}.*\Pi_{ij}
% The cvx is applied to solve the convex optimization problem above.

%% Define the two probability distributions p_{i} and q_{j} and cooridinates X(p), Y(q)
% Consider p_{i} and q_{j} are two discretes distributions 
% p_{i}: i = 1, 2,..., n1
% q_{j}: j = 1, 2,..., n2

clear all
clc

% Useful functions 
normalize = @(a)a/sum(a(:)); % Normalization function 
%% Option 2: Gaussian distributions p,q in histogram
n1 = 40;  % (not exceed 1000)
n2 = 60;  % (not exceed 1000)

x = (1:n1)'/n1; % x \in [-1,0]
y = (1:n2)'/n2;     % y \in [ 0,1]

% Useful functions 
normalize = @(a)a/sum(a(:));                           % Normalization function 
Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) ); % Gaussian distribution generator
sigma0 = 0.01; sigma1 = .02; sigma2 = .03; sigma3 = .04; sigma4 = .05;


sigma = .1;
p = Gaussian(x,.2,0.08);
q = Gaussian(y,0.8,0.1);

vmin = .02;
p = normalize( p+max(p)*vmin);
q = normalize( q+max(q)*vmin);
% p = normalize(p);
% q = normalize(q);

figure(1)
subplot(2,1,1);
bar(x, p, 'k'); axis tight;
subplot(2,1,2);
bar(y, q, 'k'); axis tight;



%% Define an linear operator for the constraints of the marginals (Coupling)
[Y,X] = meshgrid(x,y);
 C = (X+Y).^2;
% C = abs((X-Y).^(4));
% C = X.^2./Y;

figure(2)
cost = surf(X,Y,C);
cost.EdgeColor = 'none';
colorbar
%% CVX solver 
cvx_begin
cvx_solver mosek
    variable ga(n1,n2)
    minimize(sum(sum(C'.*ga)))
    subject to
     ga*ones(n2,1) == p;
     ga'*ones(n1,1) == q;
     ga >= 0;
cvx_end

% Constraint check (suppose to be 0):
sum(ga*ones(n2,1)-p)
sum(ga'*ones(n1,1)-q)

%% Map 
figure(3)
coupling = surf(X',Y',ga);
coupling.EdgeColor = 'none';
coupling.FaceColor = 'flat';
xlabel('q')
ylabel('p')
title('Coupling \Pi(p,q)')
axis tight