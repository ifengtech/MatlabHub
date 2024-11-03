% Illustrate animation and movies in Matlab.
clear
clc
close all

%% Step 1: Generate Data.
format compact;

% % spiral line
% t = linspace(0, 2*pi, 100);
% x = 5*cos(t);
% y = 2*sin(t);
% z = t;

% Lorenz system curve
t = linspace(0, 2*pi, 100);
sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
[t,a] = ode45(f,[0 100],[1 1 1]);     % Runge-Kutta 4th/5th order ODE solver
x = a(:,1);
y = a(:,2);
z = a(:,3);

%% Step 2: Draw/Render the scenario.
fig = figure;

for k=1:length(t)
    clf
    
    t_k = t(k);
    x_k = x(k);
    y_k = y(k);
    z_k = z(k);
    
    %Plot the current location of the particle.
    plot3(x_k, y_k, z_k, 'go', 'LineWidth', 3, 'MarkerSize', 10)
    
    %Plot the entire trajectory
    hold on
    plot3(x, y, z, 'b-', 'LineWidth',.5)
    
    % Decorate the plot.
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['Particle at time = ', num2str(t_k, 2), ' seconds.'])
    %view([30, 35])
    view([30+20*t_k, 35])   % auto-adjuest(rotate) the view
    
    %force Matlab to draw the image at this point.
    %drawnow
    % pause(0.2)
    videoVector(k) = getframe(fig, [10 10 520 400]);
    % videoVector(k) = getframe(fig);
    
    %% Step 4: Advance time.
    % Dont need to do anything.
    
end

%% Step 5: Save to video.
% myWriter = VideoWriter('curve');    %avi video
myWriter = VideoWriter('curve', 'MPEG-4');
myWriter.FrameRate = 20;

% Open the video object, write the video, and close the file.
open(myWriter);
writeVideo(myWriter, videoVector);
close(myWriter)
