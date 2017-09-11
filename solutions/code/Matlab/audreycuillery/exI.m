%% Exercises set I
% PhD course on Sequential Monte Carlo methods 2017

%% I.1. Importance Sampling

% (a) not a good importance sampler because the "tails" of the proposal
% distribution are not going far enough for the target distribution.
% (support of target is not included in support of proposal)

% (b) implementation of proposed importance sampler
% target : U(0,4)
% proposal : N(0,1)
% number of particules : N = 10000

N  = 10000;
s  = randn(N,1);                     % sample from proposal

% define pdf for normal distribution
pdfnorm = @(x,mu,sig2) 1/sqrt(2*pi*sig2).*exp(-0.5.*(x-mu.*ones(size(x))).^2./sig2); 
% define pdf for uniform distribution
pdfunif = @(x,min,max) (x>min & x<max)./(max-min);

wu = pdfunif(s,0,4)./pdfnorm(s,0,1); % unnormalized weights
W  = sum(wu);                        % normalizing factor
w  = wu./W;                          % normalized weights

figure(1)
clf
% stem(s,w,'Marker','none')

% kernel density estimate
[f,xi] = ksdensity(s,'Weights',w);
plot(xi,f,'k');
hold on

% plot target and proposal densities
x    = linspace(-5,5,1000);
targ = pdfunif(x,0,4);
prop = pdfnorm(x,0,1);
plot(x,targ,'b');
plot(x,prop,'r');

% weighted histogram
[H,V] = histwc(s,w,100);
bar(V,H)

legend('kernel density estimate','target','proposal')

% To improve the sampler: changing the proposal density, take a larger
% value of the variance for example, to cover the all space of the target
% density. 

%%
% (c)

Nsimu = 100;
mean_est = zeros(Nsimu,1);

for i = 1:Nsimu
   s = sqrt(2).*randn(N,1)+2;
   w = pdfunif(s,0,4)./pdfnorm(s,2,2); % no need to normalize since it's the true target density used !!! but !!! Don't forget 1/N in mean computation
%    W  = sum(wu);
%    w  = wu./W;
   
   % mean computation
   mean_est(i) = 1/N.*sum(w.*s);     % 1/N because no normalization step
end

figure(2)
clf
plot(mean_est,'o')

%%
% (d) See paper

% (e) Computation of the normalizing constant with the Monte Carlo
% approximation to see thate it is unbiased (ie always the same for each
% simuation)

Z_est = zeros(Nsimu,1);

for i = 1:Nsimu
    s  = sqrt(2).*randn(N,1)+2;
    wu = (s>0&s<4)./pdfnorm(s,2,2);
    
    % Z estimation
    Z_est(i) = 1/N*sum(wu);
end

figure(3)
clf
plot(Z_est,'o')

%%
% (f)

N = 10;

mean_est = zeros(Nsimu,1);

for i = 1:Nsimu
    
    mu   = 0;
    sig2 = 9;
    
    % samples from proposal
    s = sqrt(sig2).*randn(N,1)+mu;
    
    % compute unnormalized weights
    wu = (s>0 & s<4)./pdfnorm(s,mu,sig2);
    
    % normalizing constant estimation
    Z  = 1/N*sum(wu);
    
    % compute the mean
    mean_est(i) = 1/(N*Z).*sum(s.*wu);
    
end

figure(4)
clf
plot(mean_est,'o')

% (g) Solution to the previous corresponds to simply normalizing the sample
% weights because the 1/N simplify

%% I.2 Importance sampling in higher dimension

D_range = 1:10;
prop    = zeros(size(D_range));

N = 10000;

for D = D_range
    
    % cumpute samples from N(0,I_D) --> I_D identity matrix of dimension D
    s = randn(N,D);
    
    % compute weights
    w = (sum((s>-0.5.*ones(N,D) & s<0.5.*ones(N,D)),2)~=0)./prod(normpdf(s,0,1),2);
    
    % number of weights equal to zero
    n_w0 = sum(w==0);
    
    % proportion of weights equal to zero
    prop(D) = n_w0/N;
    
end

figure(5)
clf
prop = prop.*100;
stem(D_range,prop)
xlabel('D (dimension)')
ylabel('proportion of samples w=0')

% It seems to be an exponantial decrease of the number of samples whichs
% laid in the cube we want sample from to have an idea of the target
% density

%% I.3 An important numerical aspect

D = 1000;

% (a) 
% target density N(0,I_D)
% proposal N(0,4)

N = 10;

s = randn(N,D);

% compute target value for samples (here named p, pi in paper)
p = prod(normpdf(s,0,1),2);          

% compute proposal value for samples
q  = prod(normpdf(s,0,2),2);

% compute the weight
w = p./q;

% We obtaine NaN values for the weight because we multiply a lot of small
% values all together.

%%
% (b) Solution: to use logarithm weights
log_pi = sum(log(normpdf(s,0,1)),2);
log_q  = sum(log(normpdf(s,0,2)),2);
log_w  = log_pi-log_q;

% It's numerically better, the log-weight is computed

%%
% When \tilde{w} is too small, even the exponential of the log is still
% smaller than what the system can represent. To obtain a normalize version
% of the weights, explore the trick of computing w = exp(log(w)-max(log(w)))

C = max(log_w);
w_trick = exp(log_w - C*ones(size(log_w)));

% This approach is valid because in the normalization step, the constant
% max(log(w)) vanishes (is simplified)

%% I.4 Boostrap partcile filter for the stochastic volatility model

% dynamics model
% xt | xt-1 \sim N(phi*xt-1,sigma^2)
% observation model
% yt | xt \sim N(0,beta^2 exp(xt))

% parameter vector theta
% theta = [phi,sigma,beta]
theta = [0.98,0.16,0.70];

phi  = theta(1);
sig  = theta(2);
beta = theta(3);

T = 500 + 1; % number of time steps (+1 for initialization)
N = 500;     % number of particles

y = load('seOMXlogreturns2012to2014.csv');

%%
% reasonable assumption about the initial state x0
% --> uniform between 0 and 1

% % % Bootstrap particle filter % % %

% !!! T here is increased of one for initialization step

x = zeros(N,T);
a = zeros(N,T); % first column not used (no ancestors at initialization)
w = zeros(N,T);

% % % Initialization t = 0

% Sample from prior
x(:,1) = rand(N,1);

% Set initial weights
w(:,1) = 1/N.*ones(N,1);

% % % for t = 1 to T
for t = 2:T

    % 1. Resample and set sample ancestor indice

    % cumulative weights
    c    = zeros(N,1);
    c(1) = w(1,t-1);
    for i = 2:N
        c(i) = c(i-1) + w(i,t-1);
    end
    c = c./c(N);

    % systematic resampling
    u = 1/N*rand;
    i = 1;

    for j = 1:N
        while u>c(i)
            i = i+1;
        end
        a(j,t) = i; 
        u = u + 1/N;
    end

    % 2. Propagate
    x_prev = x(a(:,t),t-1);                   % chosen ancestors to compute the predicted ones from
    x(:,t) = sig.*randn(N,1) + phi.*x_prev;

    % 3. Compute weights
    w(:,t) = normpdf(y(t-1).*ones(N,1),zeros(N,1),beta.*sqrt(exp(x(:,t))));
    W = sum(w(:,t));     % normalization constant
    w(:,t) = w(:,t)./W;  % normalization

end

%% 
figure(6)
clf
hold on

% recover trajectories from ancestors
for n = 1:N
    
    traj     = zeros(1,T);
    traj(T)  = x(n,T);
    ancestor = a(n,T);
    for t = T-1:-1:1
        traj(t) = x(ancestor,t);
        ancestor = a(ancestor,t);
    end
        
plot(traj)

end

% mean estimation
mean_est = mean(x,1);
plot(mean_est,'r-o','MarkerFaceColor','r')
