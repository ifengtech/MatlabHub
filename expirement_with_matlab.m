% Iteration

%% 1. Golden spir
% equation: x = sqrt(1+x)
% phi: (1+sqrt(5))/2=1.618
% phi rate:(phi-1)/1=0.618

x=3;y=0;
while abs(x-y) >eps(x)
    y=x;
    x=sqrt(x+1)
end


% Plot golden solution
phi = x
x = -1:.02:4;
y1 = x; y2 = sqrt(1+x);
plot(x,y1,'-',x,y2,'-',phi,phi,'o')

%%
f = @(x) tan(sin(x)) - sin(tan(x));
ezplot(f, [-pi, pi])

%% fibnum increasing rate plot
n = 40; f = Fibnum(n);
x = 2:n;
y = f(x)./f(x-1);
plot(x,y,'-')

%% fibonacci number normal plot
n = 18;plot(1:n,fibonacci(n), '-')
% Semi-log scale( y=lnx(t) ) plot 
semilogy(fibonacci(18),'-o')


%% 
syms phi;
A = [[1,1];[phi, 1-phi]]; b= [1;1];
C = A \ b

%% Biorhythms.
t0 = datenum('1991-08-19')      % My Birthday.
t1 = fix(now);
t = (t1-28):1:(t1+28);
y = 100*[sin(2*pi*(t-t0)/23)
         sin(2*pi*(t-t0)/28)
         sin(2*pi*(t-t0)/33)];
plot(t,y)

%%
n = 2;
A = [0.99 0.01; -0.01 1.01];
A^n
[V,D]=eig(A)

%% Google PageRank
% init simple data
i = [2 6 3 4 4 5 6 1 1];
j = [1 1 2 2 3 3 3 4 6];
n = 6; G = sparse(i,j,1,n,n);
c = full(sum(G))
% 
p = 0.85; k = find(c~=0);
D = sparse(k,k, 1./c(k),n,n);
e = ones(n,1); I = speye(n,n);
x = (I - p*G*D)\e, x = x/sum(x); x', bar(x)


%% Exponential Test
t = 0:0.01:2;h = 0.00001; y=2.^t;
ydot = (2.^(t+h)-2.^t)/h; 
figure(1), hold on
subplot(1,2,1), plot(t,[y;ydot])
subplot(1,2,2), plot(t,ydot./y); axis([0 2 0.5 0.9])
hold off

% function s = expex(t)
% % calc exp on t
%     s=1;term = 1; n = 0;r=0;
%     while r ~= s
%         r = s; n=n+1;term = (t/n)*term; s = s+term;
%     end
% end


%% 2D firgure ploting with Complex Number
figure(2);
theta = (1:2:17)'*pi/8;
z = exp(i*theta),subplot(2,2,1), p1 = plot(z);
theta = (0:3:15)'*(2*pi/5)+pi/2;
z = exp(i*theta),subplot(2,2,2), p2 = plot(z);

set(p,'linewidth', 4, 'color', 'blue');
set(p,'linewidth', 4, 'color', 'red');
axis square, axis off


%% Magic square
% case 1: odd n.
n = 5;
[I,J] = ndgrid(1:n);
A = mod(I+J+(n-3)/2,n);B = mod(I+2*J-2,n);
M = n*A + B +1; disp(M)

% case 2: doubly-even n.
n = 8;
M = reshape(1:n^2,n,n)';
[I,J] = ndgrid(1:n);
K = fix(mod(I,4)/2) == fix(mod(J,4)/2);
M(K) = n^2 +1 -M(K); disp(M)

%% Lorenz System Plotting
% Solve over time interval [0,100] with initial conditions [1,1,1]
% ''f'' is set of differential equations
% ''a'' is array containing x, y, and z variables
% ''t'' is time variable

sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
[t,a] = ode45(f,[0 100],[1 1 1]);     % Runge-Kutta 4th/5th order ODE solver
plot3(a(:,1),a(:,2),a(:,3))
