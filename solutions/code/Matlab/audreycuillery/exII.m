%% Exercises set II
% PhD course on Sequential Monte Carlo methods 2017

%% II.1 Likelihood estimates for the stochastic volatility model

% (b) likelihood estimation depending on beta

% dynamics model
% xt | xt-1 \sim N(phi*xt-1,sigma^2)
% observation model
% yt | xt \sim N(0,beta^2 exp(xt))

% parameter vector theta
% theta = [phi,sigma,beta]
theta = [0.98,0.16,0.70];    % from previous exercise I.4

phi  = theta(1);
sig  = theta(2);
% !!! beta has to be estimated

T = 500 + 1; % number of time steps (+1 for initialization)
N = 500;     % number of particles

y = load('seOMXlogreturns2012to2014.csv');

%%

beta_grid = 0.4:0.2:2;          % grid to test beta values
n_beta    = length(beta_grid);  % number of beta tested
n_simu    = 10;                 % number of simulations to run for each value of beta  

lik_est   = zeros(n_simu,n_beta);

% % % % % Likelihood estimation % % % % %

for ibeta = 1:n_beta

    beta = beta_grid(ibeta);

    for isimu = 1:n_simu

        % % % Bootstrap particle filter % % %

        % !!! T here is increased of one for initialization step
        
        x  = zeros(N,T);
        a  = zeros(N,T); % first column not used (no ancestors at initialization)
        wu = zeros(N,T); % unnormalized weights (to compute estimate of the likelihood)
        w  = zeros(N,T);
        
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
            wu(:,t) = normpdf(y(t-1).*ones(N,1),zeros(N,1),beta.*sqrt(exp(x(:,t)))); % previous weight does not appear because of resampling step which sets weight to 1/N
            W = sum(wu(:,t));     % normalization constant
            w(:,t) = w(:,t)./W;   % normalization

        end % end loop over time for boostrap filter
        
        lik_est(isimu,ibeta) = sum(log(sum(wu(:,2:end),1))-log(N).*ones(1,T-1),2);

    end % end loop over simu

end % end loop over beta

%% 
figure(7)
clf
hold on

boxplot(lik_est,beta_grid)
xlabel('beta parameter value')
ylabel('estimated log likelihood')

%%
% (b) influence of N and T over the log-likelihood esttimate

% INFLUENCE OF N

N_grid = [10,100,500,1000,2000];
n_N    = length(N_grid);

lik_est   = zeros(n_simu,n_beta,n_N);

for iN = 1:n_N
    
    N = N_grid(iN);

    % % % % % Likelihood estimation % % % % %

    for ibeta = 1:n_beta

        beta = beta_grid(ibeta);

        for isimu = 1:n_simu

            % % % Bootstrap particle filter % % %

            % !!! T here is increased of one for initialization step
            
            x  = zeros(N,T);
            a  = zeros(N,T); % first column not used (no ancestors at initialization)
            wu = zeros(N,T); % unnormalized weights (to compute estimate of the likelihood)
            w  = zeros(N,T);

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
                wu(:,t) = normpdf(y(t-1).*ones(N,1),zeros(N,1),beta.*sqrt(exp(x(:,t)))); % previous weight does not appear because of resampling step which sets weight to 1/N
                W = sum(wu(:,t));     % normalization constant
                w(:,t) = w(:,t)./W;   % normalization

            end % end loop over time for boostrap filter

            lik_est(isimu,ibeta,iN) = sum(log(sum(wu(:,2:end),1))-log(N).*ones(1,T-1),2);

        end % end loop over simu

    end % end loop over beta

end % end loop over N

%%
figure(8)
clf
hold on
V = var(lik_est,1);
V = squeeze(V);
h = cell(1,n_N);    % for legend
for iN = 1:n_N
    col = rand(1,3);
    h{iN} = plot(N_grid,V(iN,:),'+','MarkerEdgeColor',col);
end
xlabel('beta parameter value')
ylabel('variance of estimated log likelihood')  
legend([h{:}],num2str(N_grid'))

%%
% INFLUENCE OF T

N = 500;

T_grid = [10,50,100,200,500];
n_T    = length(T_grid);
T_grid = T_grid + ones(1,n_T);

lik_est   = zeros(n_simu,n_beta,n_T);

for iT = 1:n_T
    
    T = T_grid(iT);

    % % % % % Likelihood estimation % % % % %

    for ibeta = 1:n_beta

        beta = beta_grid(ibeta);

        for isimu = 1:n_simu

            % % % Bootstrap particle filter % % %

            x  = zeros(N,T);
            a  = zeros(N,T); % first column not used (no ancestors at initialization)
            wu = zeros(N,T); % unnormalized weights (to compute estimate of the likelihood)
            w  = zeros(N,T);

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
                x_prev = x(a(:,t),t-1);                   % chosen ancestors to compute from
                x(:,t) = sig.*randn(N,1) + phi.*x_prev;

                % 3. Compute weights
                wu(:,t) = normpdf(y(t-1).*ones(N,1),zeros(N,1),beta.*sqrt(exp(x(:,t)))); % previous weight does not appear because of resampling step which sets weight to 1/N
                W = sum(wu(:,t));     % normalization constant
                w(:,t) = w(:,t)./W;   % normalization

            end % end loop over time for boostrap filter

            lik_est(isimu,ibeta,iT) = sum(log(sum(wu(:,2:end),1))-log(N).*ones(1,T-1),2);

        end % end loop over simu

    end % end loop over beta

end % end loop over T

%%
figure(9)
clf
hold on
V = var(lik_est,1);
V = squeeze(V);
h = cell(1,n_T);    % for legend
for iT = 1:n_T
    col = rand(1,3);
    h{iT} = plot(T_grid,V(iT,:),'+','MarkerEdgeColor',col);
end
xlabel('beta parameter value')
ylabel('variance of estimated log likelihood')  
legend([h{:}],num2str(T_grid'))

%% II.2 Fully adapted particle filter

% the locally optimal proposals (p(xt|xt-1,yt) and p(yt|xt-1)) can be
% computed when Yt|Xt is conjugate to Xt|Xt-1

% (a) for each following model, possible to implement fully adapted
% particle filter ?

% (i)   Yes 
% (ii)  Yes
% (iii) No

% (b) Fully adapted particle filter for model (ii)

% xt+1 = cos(xt)^2 + vt        vt \sim N(0,1)      P = 1;
%  yt  = 2*xt      + et        et \sim N(0,0.01)   R = 0.01;

% reformulation for notation:
% xt+1 = m(xt) + vt        vt \sim N(0,1)      P = 1;
%  yt  = H*xt  + et        et \sim N(0,0.01)   R = 0.01;

m = @(x) cos(x).^2;
P = 1;
H = 2;
R = 0.01;

T = 100;

% % % data creation
% true trajectory and observations
x0_tr   = 0;
x_tr    = zeros(1,T);
x_tr(1) = x0_tr;

y    = zeros(1,T);
y(1) = H*x0_tr + sqrt(R)*randn;

for t = 2:T

    x_tr(t) = m(x_tr(t-1)) + sqrt(P)*randn;
    y(t)    = H*x_tr(t)    + sqrt(R)*randn;
    
end

figure(10)
clf
plot(1:T,x_tr)
title('true trajectory')
xlabel('time')
ylabel('x')

%%
N = 500; % number of particles
T = 100;
T = T+1; 
% !!! T here is increased of one for initialization step

% % % Fully adapted particle filter % % %

x = zeros(N,T);
a = zeros(N,T); % first column not used (no ancestors at initialization)

% NB: no need weights vector because 1/N each time

% what we have:
%    - p(xt|xt-1) = N(xt|cos(xt-1)^2,P) = N(xt|m,P)     (notation)
%    - p(yt|xt)   = N(yt|2*xt,R)        = N(yt|H*xt,R)  (notation)
% what we want: 
%    - p(xt|xt-1,yt)
%    - p(yt|xt-1)
% we know:
%              p(xt|xt-1)p(yt|xt) =  p(xt|xt-1,yt)p(yt|xt-1)
% conjugaison lemma between normal distributions gives: 
%    - p(xt|xt-1,yt) = N(xt|m'(xt-1,yt),P')
%    - p(yt|xt-1)    = N(yt|Hm,HPHt+R)
% with:
%    - H = 2                         (see observation model)
%    - m = cos(xt-1)^2               (see observation model)
%    - K = PHt(HPHt+R)^(-1)          --> is a constant ! 
%    - m'(xt-1,yt) = m + K(yt-H*m) 
%    - P' = (I-KH)P                  --> is a constant !                 

K  = H*P/(H^2*P+R);
Pp = (1-H*K)*P;

% define pdf resample = p(yt|xt-1) = N(yt|Hm,HPHt+R)
sig2 = H^2*P+R;
pdfresample = @(y,x) 1/sqrt(2*pi*sig2).*exp(-0.5*(y-H*m(x)).^2./sig2);

% % % Initialization t = 0

% Sample from prior
x(:,1) = randn(N,1);

% % % for t = 1 to T
for t = 2:T

    % 1. Resample and set sample ancestor indice

    w = pdfresample(y(t-1),x(:,t-1));
    
    % cumulative weights
    c    = zeros(N,1);
    c(1) = w(1);
    for i = 2:N
        c(i) = c(i-1) + w(i);
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
    x_prev = x(a(:,t),t-1);                      % chosen ancestors to compute the predicted ones from
    % m'(xt-1,yt) = m + K(yt-H*m)
    mp     = m(x_prev) + K.*(y(t-1)-H*m(x_prev));
    x(:,t) = sqrt(Pp).*randn(N,1) + mp;

    % 3. Compute weights
    % no need because 1/N for each particle

end

% display
figure(11)
clf
hold on
title('fully adapted particle filter')
xlabel('time')
ylabel('x')

plot(2:T,x_tr,'b.')      % true trajectory
plot(1:T,mean(x,1),'ko') % mean of particle filter

for n = 1:N
    
    traj     = zeros(1,T);
    traj(T)  = x(n,T);
    ancestor = a(n,T);
    for t = T-1:-1:1
        traj(t) = x(ancestor,t);
        ancestor = a(ancestor,t);
    end
        
plot(traj,'r')

end

plot(2:T,x_tr,'b.')      % true trajectory
plot(1:T,mean(x,1),'ko') % mean of particle filter

legend('true trajectory','mean of fully adapted PF','particle trajectories')

%%
% % % Bootstrap particle filter % % %

likelihood = @(y,x) 1/sqrt(2*pi*R).*exp(-0.5.*(y-H.*x).^2./R);

% !!! T here is increased of one for initialization step

x_boot = zeros(N,T);
a_boot = zeros(N,T); % first column not used (no ancestors at initialization)
w_boot = zeros(N,T);

% % % Initialization t = 0

% Sample from prior
x_boot(:,1) = randn(N,1);

% Set initial weights
w_boot(:,1) = 1/N.*ones(N,1);

% % % for t = 1 to T
for t = 2:T

    % 1. Resample and set sample ancestor indice

    % cumulative weights
    c    = zeros(N,1);
    c(1) = w_boot(1,t-1);
    for i = 2:N
        c(i) = c(i-1) + w_boot(i,t-1);
    end
    c = c./c(N);

    % systematic resampling
    u = 1/N*rand;
    i = 1;

    for j = 1:N
        while u>c(i)
            i = i+1;
        end
        a_boot(j,t) = i; 
        u = u + 1/N;
    end

    % 2. Propagate
    x_prev = x_boot(a_boot(:,t),t-1);                   % chosen ancestors to compute the predicted ones from
    x_boot(:,t) = sqrt(P).*randn(N,1) + m(x_prev);

    % 3. Compute weights
    w_boot(:,t) = likelihood(y(t-1),x_prev);
    W = sum(w_boot(:,t));          % normalization constant
    w_boot(:,t) = w_boot(:,t)./W;  % normalization

end

% display
figure(12)
clf
hold on
title('bootstrap particle filter')
xlabel('time')
ylabel('x')

plot(2:T,x_tr,'b.')                % true trajectory
plot(1:T,sum(w_boot.*x_boot),'ko') % mean of particle filter

for n = 1:N
    
    traj     = zeros(1,T);
    traj(T)  = x_boot(n,T);
    ancestor = a_boot(n,T);
    for t = T-1:-1:1
        traj(t) = x_boot(ancestor,t);
        ancestor = a_boot(ancestor,t);
    end
        
plot(traj,'r')

end

plot(2:T,x_tr,'b.')                % true trajectory
plot(1:T,sum(w_boot.*x_boot),'ko') % mean of particle filter

legend('true trajectory','mean of fully adapted PF','particle trajectories')

%%
% variance comparison between boostrap filter and fully adapted particle
% filter

m_fapf   = mean(x);
var_fapf = var(x);
m_boot   = sum(w_boot.*x_boot);
var_boot = zeros(1,T);
for i = 1:T
var_boot(i) = var(x_boot(:,i),w_boot(:,i));
end

nb_sig = 1; % number of sigma to plot

figure(13)
clf
hold on

upp_fapf = m_fapf + nb_sig.*sqrt(var_fapf);
low_fapf = m_fapf - nb_sig.*sqrt(var_fapf);

upp_boot = m_boot + nb_sig.*sqrt(var_boot);
low_boot = m_boot - nb_sig.*sqrt(var_boot);

filled  = [upp_boot,fliplr(low_boot)];
xpoints = [1:T,fliplr(1:T)];
fillhandle = fill(xpoints,filled,'g');
set(fillhandle,'EdgeColor','g','FaceAlpha',0.4,'EdgeAlpha',1);

filled  = [upp_fapf,fliplr(low_fapf)];
xpoints = [1:T,fliplr(1:T)];
fillhandle = fill(xpoints,filled,'c');
set(fillhandle,'EdgeColor','c','FaceAlpha',0.4,'EdgeAlpha',1);

plot(2:T,x_tr,'b-.')                % true trajectory

nb_sig_str = num2str(nb_sig);
title(strcat('comparison with ', nb_sig_str ,' standard deviation'))
xlabel('time')
ylabel('x')

%% II.3 Likelihood estimator for the APF

% see paper

%% II.4 Forgetting

% boostrap particle filter for LGSS model from H.2 (a)-(c) but with Q=0



