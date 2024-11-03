%
clear;
clc;
close all;

%% Plotting One-Dot of variable solution.
t0 = 0;
tfinal = 10;
X0 = 1;

[t, dotx] = ode45(@(t,x) 2*x, [t0 tfinal], X0);

plot(t,dotx)
grid


%% Plotting One-Dot of vector solution.
clear;
clc;
close all;

t0 = 0;
tfinal = 10;
X0 = [0; 1; 3; 3; 5; 6;];

D = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;];
A = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0;];
L = D - A;  % coeffect

[t, dotX] = ode45(@(t,X) -L*X, [t0 tfinal], X0);

plot(t,dotX(:,1), t,dotX(:,2), t,dotX(:,3), t,dotX(:,4), t,dotX(:,5), t,dotX(:,6))
grid

%% Plotting Two-Dot of variable solution.
% convert two-dot equation into two demensions vector of one-dot equation.
[t,y] = ode45(@(t,y) [y(2); (1-y(1)^2)*y(2)-y(1)],[0 20],[2; 0]);

% plot $y_1$ and $y_2$ curves.
plot(t,y(:,1),'-o',t,y(:,2),'-o')
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')
grid on



