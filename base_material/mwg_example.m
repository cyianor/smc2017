function gibbs_nonlinear
close all
%% Plot target
xv = -1.5:0.01:1.5;
yv = -2:0.01:2;
[xx,yy] = meshgrid(xv,yv);
a = 1;
b = 0.1;
z = logpi([xx(:) yy(:)]', a,b);
z = reshape(z, length(yv), length(xv));
contour(xv,yv,exp(z));
axis xy

%% Run sampler
M = 150;
X = zeros(1,M);
Y = zeros(1,M);
sq = 0.5;
hold on;
plot(X(1),Y(1),'r.','MarkerSize',12);
xlabel('x_1');
ylabel('x_2');
set(gca,'FontSize',12);
pause;
ar = 0;
for(m = 2:M)
    Xp = X(m-1) + sq*randn(1);
    ap = logpi([Xp;Y(m-1)],a,b) - logpi([X(m-1);Y(m-1)],a,b);
    if(rand(1) < exp(ap))
        X(m) = Xp;
        ar = ar+1;
        
        plot([X(m-1) X(m)],[Y(m-1) Y(m-1)],'r.','MarkerSize',12);
    else
        X(m) = X(m-1);
        plot(Xp,Y(m-1),'kx','linewidth',1,'Markersize',8);      
    end
    pause(0.1)
    
    Y(m) = X(m)^3 + sqrt(b)*randn();
    plot([X(m) X(m)],[Y(m-1) Y(m)],'r.','MarkerSize',9);
    pause(0.1)
end

end
%----------------------------------
function z = logpi(x,a,b)
z = 1/a*(2*x(1,:)+sin(x(1,:)*2*pi)).^2 + 1./(b) .* (x(2,:) - x(1,:).^3).^2;
z = -1/2*z;
end