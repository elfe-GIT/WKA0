function postprocess(sys,t,y)
% do plots and animations
%   start with some data preparation
u1X=y(:, 1);
u1Y=y(:, 2);
u1Z=y(:, 3);
u2X=y(:, 4);
u2Y=y(:, 5);
u2Z=y(:, 6);
u3X=y(:, 7);
u3Y=y(:, 8);
u3Z=y(:, 9);
v1X=y(:,10);
v1Y=y(:,11);
v1Z=y(:,12);
v2X=y(:,13);
v2Y=y(:,14);
v2Z=y(:,15);
v3X=y(:,16);
v3Y=y(:,17);
v3Z=y(:,18);

xRange = [min([u1X;u2X;u3X]),max([u1X;u2X;u3X])];
yRange = [min([u1Y;u2Y;u3Y]),max([u1Y;u2Y;u3Y])];
zRange = [min([u1Z;u2Z;u3Z]),max([u1Z;u2Z;u3Z])];

% body-radius
R = (sys.theta/max(sys.theta)).^(1/3);

%%
% plot trajectories
figure
p=plot3(u1X,u1Y,u1Z,'-ko', ...
        u2X,u2Y,u2Z,'-ko', ...
        u3X,u3Y,u3Z,'-ko');
grid on;
p(1).MarkerSize = 10*R(1);
p(1).Color      = '#ff0000';
p(2).MarkerSize = 10*R(2);
p(2).Color      = '#00ff00';
p(3).MarkerSize = 10*R(3);
p(3).Color      = '#0000ff';
title('trajectories of bodies')
xlabel('x/L ->');
ylabel('y/L ->');
zlabel('z/L ->');
legend('body 1','body 2', 'body 3');


%% 
% check total momentum (must remain constant)
I = [sys.theta(1)*v1X + sys.theta(2)*v2X + sys.theta(3)*v3X, ...
     sys.theta(1)*v1Y + sys.theta(2)*v2Y + sys.theta(3)*v3Y, ...
     sys.theta(1)*v1Z + sys.theta(2)*v2Z + sys.theta(3)*v3Z];

figure
q=plot(t,I(:,1),t,I(:,2),t,I(:,3));
title('total momentum in three spatial directions ')
xlabel('t/T ->');
ylabel('I/I[0] ->');
legend('x-direction','y-direction', 'z-direction');

%% animated plot
% generate avi ... or anything else

fid = figure;
% The initial plot (void).
plot3(0, 0, 0, '+');
axis([xRange(1) xRange(2) yRange(1) yRange(2) zRange(1) zRange(2)]);
grid on;

%% Set up the movie.
% from: https://kawahara.ca/matlab-make-a-movievideo/
writerObj = VideoWriter('T3BP.mp4', 'MPEG-4'); % Name it.
writerObj.FrameRate = 15; % How many frames per second.
open(writerObj); 

for T=1:size(t)
    figure(fid); % Makes sure you use your desired frame.
    % plot for time index T
    p=plot3( u1X(T),u1Y(T),u1Z(T),'-ko', ...
             u2X(T),u2Y(T),u2Z(T),'-ko', ...
             u3X(T),u3Y(T),u3Z(T),'-ko', ...
            [u1X(T),u2X(T)],[u1Y(T),u2Y(T)],[u1Z(T),u2Z(T)],'--', ...
            [u1X(T),u3X(T)],[u1Y(T),u3Y(T)],[u1Z(T),u3Z(T)],'--', ...
            [u2X(T),u3X(T)],[u2Y(T),u3Y(T)],[u2Z(T),u3Z(T)],'--'); 
    grid on;
    axis([xRange(1) xRange(2) yRange(1) yRange(2) zRange(1) zRange(2)]);
    % bodies
    p(1).MarkerSize = 25*R(1);
    p(1).Color      = '#ff0000';
    p(1).MarkerFaceColor = '#ff0000';
    p(2).MarkerSize = 25*R(2);
    p(2).Color      = '#00ff00';
    p(2).MarkerFaceColor = '#00ff00';
    p(3).MarkerSize = 25*R(3);
    p(3).Color      = '#0000ff';
    p(3).MarkerFaceColor = '#0000ff';
    % connecting lines for a better illustration of geometry
    p(4).Color      = '#333333';
    p(5).Color      = '#333333';
    p(6).Color      = '#333333';
    % the usual stuff ...
    xlabel('x/L ->');
    ylabel('y/L ->');
    zlabel('z/L ->');
    legend('body 1','body 2', 'body 3');
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end
close(writerObj); % saves the movie.

end

