clear all;
close all;

N = 10000;

% random number seed
sd=1;
rng(sd);

% proposal
mu = 0;
sigma = sqrt(1);


% target
a = 0;
b = 4;
x = a + (b-a).*rand(N,1);

%y = norminv(x,mu,sigma);
y = sigma.*(-sqrt(2)*erfcinv(2*x))+mu;

figure(1);
histogram(x);

figure(2)
histogram(y,'Normalization','pdf');
hold on
g = -5:0.1:5;
f = exp(-(g-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(g,f,'LineWidth',2);

% Estimation of variance

%var = 1/N*sum((y-mu).^2)
