function postprocess(sys,t,y,range)
% do plots
%   start with some data preparation

% tEnd = 0 means: plot it all
if range(2)==0
   range(2)=t(end); 
end

% retrieve indeces of start/end-values
% tower
QT = y(:, 1: 3);
PT = y(:,15:17);
% shaft
QS = y(:, 4: 5);
PS = y(:,18:19);
% blades
QB = y(:, 6:14);
PB = y(:,20:28);

r = [find(abs(t-range(1))==min(abs(t-range(1)))),...
     find(abs(t-range(2))==min(abs(t-range(2))))];


%trim results space
t  = t(r(1):r(2));
QT = QT(r(1):r(2), :);
QS = QS(r(1):r(2), :);
QB = QB(r(1):r(2), :);
PT = PT(r(1):r(2), :);
PS = PS(r(1):r(2), :);
PB = PB(r(1):r(2), :);

%+++++++++++++++++++++++++++++++++++++++++++++++++
subplot(2,2,1);

%%
% plot tower
plot(t,QT);
title('tower');
xlabel('t/s ->');
ylabel('W,Θ ->');
legend('W_1','W_2','Θ');

subplot(2,2,2);
%%
% plot shaft
plot(t,QS);
title('shaft')
xlabel('t/s ->');
ylabel('V ->');
legend('V_1','V_2');

subplot(2,2,3);
%%
% plot relative velocities
plot(t,QB)
title('blades')
xlabel('t/s ->');
ylabel('u ->');
legend('u_{11}','u_{12}','u_{13}','u_{21}','u_{22}','u_{23}','u_{31}','u_{32}','u_{33}');

subplot(2,2,4);
%%
% plot contact forces
plot(QB(:,:),PB(:,:));
title('phase portrait blades')
xlabel('QB ->');
ylabel('PB ->');
legend('u_{11}','u_{12}','u_{13}','u_{21}','u_{22}','u_{23}','u_{31}','u_{32}','u_{33}');
%ylim([-5 25])
%set(gcf,'Color','white')
end

