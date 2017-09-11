%% Exercises set III
% PhD course on Sequential Monte Carlo methods 2017

%% III.1 Metropolis-Hastings

% samples from pi(x) \propto sin2(x)exp(-|x|) x in R
% gaussian random walk as proposal q(x|x') = N(x|x',sig2)
% plot results histogram with pi overlaid
% different values of sig2

p = @(x) sin(x).^2.*exp(-(abs(x))); % proportional value of pi

X = linspace(-10,10,1000);
Y = p(X);
Y = Y./max(Y);

figure(1)
clf
plot(X,Y)
hold on

M    = 1000;      % number of MCMC steps
sig2 = 0.5;          % variance proposal
q    = @(x,mu) 1/(sqrt(2*pi*sig2))*exp(-0.5*(x-mu)^2/sig2);
x    = 1;          % initialization
s    = zeros(M,1); % stored samples
s(1) = x;

for i = 1:M
   
    x_star = sqrt(sig2)*randn + x;
    ratio  = p(x_star)*q(x,x_star)/(p(x)*q(x_star,x));
    alpha  = min(1,ratio);
    
    u = rand;
    if u <= alpha
        x = x_star;
    end
    
    s(i) = x;
    
end

nbins = 100;
[H,C] = hist(s,nbins);
H = H./max(H);
bar(C,H)

%%
% test with different values of sig2

figure(2)
clf 
hold on
ax1 = subplot(2,2,1);
bar(ax1,C,H)
hold on
plot(X,Y)
title('sig2 = 1')

sig2_grid = [0.1 2 5];          % variance proposal

for j = 1:3
    
    sig2 = sig2_grid(j);
    q    = @(x,mu) 1/(sqrt(2*pi*sig2))*exp(-0.5*(x-mu)^2/sig2);
    
    x    = 1;          % initialization
    s    = zeros(M,1); % stored samples
    s(1) = x;
    
    for i = 1:M

        x_star = sqrt(sig2)*randn + x;
        ratio  = p(x_star)*q(x,x_star)/(p(x)*q(x_star,x));
        alpha  = min(1,ratio);

        u = rand;
        if u <= alpha
            x = x_star;
        end

        s(i) = x;

    end

    [H,C] = hist(s,nbins);
    H = H./max(H);
    
    ax1 = subplot(2,2,j+1);
    bar(ax1,C,H)
    hold on
plot(X,Y)
    title(strcat('sig2 =',num2str(sig2)))

end

%% III.2 Gibbs sampling

% pi(x) = N([7;3],[0.3 0.1;... 
%                  0.1  1])

M    = 1000;

mu   = [7;3];
cov  = [0.3 0.1;...
       0.1  1 ];

x0     = [0;0];
x      = zeros(2,M);
x(:,1) = x0;

sig2_x1x2 = cov(1,1)-(cov(1,2)^2/cov(2,2));
sig2_x2x1 = cov(2,2)-(cov(2,1)^2/cov(1,1));

for i = 2:M
    
    mu_x1x2 = mu(1) + cov(1,2)/cov(2,2)*(x(2,i-1) - mu(2));
    x(1,i)  = sqrt(sig2_x1x2)*randn + mu_x1x2;
    mu_x2x1 = mu(2) + cov(2,1)/cov(1,1)*(x(1,i)   - mu(1));
    x(2,i)  = sqrt(sig2_x2x1)*randn + mu_x2x1;

end

figure(3)
clf
plot(x(1,1:10),x(2,1:10),'-')
hold on
plot(x(1,:),x(2,:),'.')

C     = chol(cov);
angle = linspace(0,2*pi,200)';
xy    = [cos(angle) sin(angle)];
XY    = xy*C;
plot(mu(1)+XY(:,1), mu(2)+XY(:,2), 'r-', ...
     mu(1)+2*XY(:,1), mu(2)+2*XY(:,2), 'r--')

%% III.3 Resampling

% sample from a distribution of our choice, and provide (positive) weights
% and normalize them
% estimate of mean (with sum of weighted samples)
% estimate of mean after some resampling schemes: multinomial, stratified,
% systematic
% even though the resampling is unbiased, the variance of these estimators
% is always larger than (or possibly equal to) the variance of the
% estimator of weighted samples --> resampling adds variance

p = @(x) 0.3*1/sqrt(2*pi*0.1).*exp(-0.5.*(x-2).^2./0.1) + 0.7*1/sqrt(2*pi*1).*exp(-0.5.*(x-0.5).^2./2);

mu_q   = 0.5;
sig2_q = 4;
q      = @(x) 1/sqrt(2*pi*sig2_q).*exp(-0.5.*(x-mu_q).^2./sig2_q);

x = linspace(-10,10);
y = p(x);

figure(4)
clf
plot(x,y)
hold on
y = q(x);
plot(x,y,'r')

%%
% One simulation (to plot resampling schemes)

N = 100; % number of samples

% sample from proposal
s = sqrt(sig2_q).*randn(N,1) + mu_q;
% weight
w = p(s)./q(s);
% normalize weights
W = sum(w);
w = w./W;

% % weighted histogram
% [H,V] = histwc(s,w,100);
% bar(V,H)

% estimate mean 
m = sum(w.*s);

% cumulative sum weights
c    = zeros(N,1);
c(1) = w(1);

for i = 2:N
    c(i) = c(i-1)+w(i);
end

% multinomial resampling
ind = zeros(N,1);
wm  = 1/N.*ones(N,1);

figure(4)
clf
hold on
title('multinomial resampling')
for i = 1:N
    u = rand;
    ind(i) = find(u*ones(N,1)<c,1);
    plot(cos(u*2*pi),sin(u*2*pi),'.')
end

sm = s(ind);

mm = 1/N*sum(sm);

% stratified resampling
ind = zeros(N,1);
wst = 1/N.*ones(N,1);

figure(5)
clf
hold on
title('stratified resampling')
for i = 1:N
    u = rand/N + (i-1)/N;
    U = u.*ones(N,1);
    ind(i) = find(U<c,1);
    plot(cos(u*2*pi),sin(u*2*pi),'.')
end

sst = s(ind);

mst = 1/N*sum(sst);

% systematic resampling
ssy = zeros(N,1);      % samples from systematic resampling
wsy = 1/N.*ones(N,1);

u = 1/N*rand;
i = 1;

figure(6)
clf
hold on
title('systematic resampling')

for j = 1:N
    while u>c(i)
        i = i+1;
    end
    ssy(j) = s(i); 
    u = u + 1/N;
    plot(cos(u*2*pi),sin(u*2*pi),'.')
end

msy = 1/N.*sum(ssy);

%%
% Several simulations to analyse variance for each resampling scheme

N     = 100;            % number of samples
Nsimu = 100;            % number of simulations

m     = zeros(Nsimu,1); % mean without resampling
mm    = zeros(Nsimu,1); % mean with multinomial resampling
mst   = zeros(Nsimu,1); % mean with stratified resampling
msy   = zeros(Nsimu,1); % mean with systematic resampling

for i_simu = 1:Nsimu

    % sample from proposal
    s = sqrt(sig2_q).*randn(N,1) + mu_q;
    % weight
    w = p(s)./q(s);
    % normalize weights
    W = sum(w);
    w = w./W;

    % % weighted histogram
    % [H,V] = histwc(s,w,100);
    % bar(V,H)

    % estimate mean 
    m(i_simu) = sum(w.*s);

    % cumulative sum weights
    c    = zeros(N,1);
    c(1) = w(1);

    for i = 2:N
        c(i) = c(i-1)+w(i);
    end

    % multinomial resampling
    ind = zeros(N,1);
    wm  = 1/N.*ones(N,1);

    for i = 1:N
        u = rand;
        ind(i) = find(u*ones(N,1)<c,1);
    end

    sm = s(ind);

    mm(i_simu) = 1/N*sum(sm);

    % stratified resampling
    ind = zeros(N,1);
    wst = 1/N.*ones(N,1);

    for i = 1:N
        u = rand/N + (i-1)/N;
        U = u.*ones(N,1);
        ind(i) = find(U<c,1);
    end

    sst = s(ind);

    mst(i_simu) = 1/N*sum(sst);

    % systematic resampling
    ssy = zeros(N,1);      % samples from systematic resampling
    wsy = 1/N.*ones(N,1);

    u = 1/N*rand;
    i = 1;

    for j = 1:N
        while u>c(i)
            i = i+1;
        end
        ssy(j) = s(i); 
        u = u + 1/N;
    end

    msy(i_simu) = 1/N.*sum(ssy);

end % end loop simulation (Nsimu)

%%
% boxplot of Nsimu simulations

figure(7)
clf
labels = {'none','multinomial','stratified','systematic'};
data   = [m,mm,mst,msy];
boxplot(data,'Labels',labels)

%% III.4 Path-space view

