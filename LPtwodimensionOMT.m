%% Multi-Marginal Optimal mass transport problem
% This script numerically solves an OMT problem as follows:
% 1. Transportation cost: C_{ijk} = \frac{x_{i}^2}{t_{k}} + \frac{y_{i}^2}{T - t_{k}}
% 2. Transportation map: Pi_{ijk} with its three marginials: p_i, q_i, r_k
%    such that \sum_{jk} = p_{i}, \sum_{ik} = q_{j}, \sum_{ij} = r_k.
%    p_i and q_j are two probability vector   
% 3. The mariginal along the t cooridinates has an upper-bound: \frac{1.1}{T}.

% The optimizer is the vector r_k which minimize the transportation cost:
% min_{ijk} C_{ijk}.*\Pi_{ijk}
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
n1 = 100;
n2 = 100;
x = (0:n1-1)'/n1;
y = (0:n2-1)'/n2+2;

Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) ); % Gaussian distribution generator

sigma = .1;
p = Gaussian(x,.2,0.14);
q = Gaussian(y,2.5,0.1);

vmin = .02;
% p = normalize( p+max(p)*vmin);
% q = normalize( q+max(q)*vmin);

p = normalize(p);
q = normalize(q);

figure(1)
subplot(2,1,1);
bar(x, p, 'k'); axis tight;
subplot(2,1,2);
bar(y, q, 'k'); axis tight;



%% Define an linear operator for the constraints of the marginals (Coupling)
[Y,X] = meshgrid(x,y);
% C = (X-Y).^2;
C = (X-Y).^(1.1);

figure(2)
cost = surf(X,Y,C);
cost.EdgeColor = 'none';
colorbar


%% CVX solver 
cvx_begin
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
map = surf(X',Y',ga);
map.EdgeColor = 'none';
xlabel('q')
ylabel('p')
title('Coupling \Pi(p,q)')
axis tight
