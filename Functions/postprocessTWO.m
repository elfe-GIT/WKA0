% postprocessTWO
% In plots with u: Lines can be toggled (visible/invisible) by clicking on the legend entries

function postprocessTWO(sys, t, y, range)

% tEnd = 0 means: plot it all
if range(2)==0
   range(2)=t(end); 
end

% retrieve indeces of start/end-values
QT = y(:, 1: 3);     % tower position
PT = y(:,15:17);     % tower momentum
QS = y(:, 4: 5);     % shaft position
PS = y(:,18:19);     % shaft momentum
QB = y(:, 6:14);     % blades position
PB = y(:,20:28);     % blades momentum

r = [find(abs(t-range(1))==min(abs(t-range(1)))), ...
     find(abs(t-range(2))==min(abs(t-range(2))))];

% Trim results
t  = t(r(1):r(2));
QT = QT(r(1):r(2), :);
QS = QS(r(1):r(2), :);
QB = QB(r(1):r(2), :);
PT = PT(r(1):r(2), :);
PS = PS(r(1):r(2), :);
PB = PB(r(1):r(2), :);

%-------------------------------------------------
fig = figure;
fig.WindowState = 'maximized';
%-------------------------------------------------
subplot(2,2,1);
plot(t,QT);
title('tower');
xlabel('t/s ->');
ylabel('W/m, Θ/rad ->');
legend('W_1','W_2','Θ','Location','best');
grid on;
%-------------------------------------------------
subplot(2,2,2);
%%
% plot shaft
plot(t,QS);
title('shaft')
xlabel('t/s ->');
ylabel('V/m ->');
legend('V_2','V_3');
%-------------------------------------------------
subplot(2,2,3);
hQB = plot(t, QB);
title('blades');
xlabel('t/s ->');
ylabel('u/m ->');
leg = legend(hQB, ...
    {'u_{11}','u_{12}','u_{13}','u_{21}','u_{22}','u_{23}','u_{31}','u_{32}','u_{33}'}, ...
    'Location','best');
grid on;

% Interaktive Legende aktivieren
set(leg, 'ItemHitFcn', @(src, event) ...
    set(event.Peer, 'Visible', toggle(event.Peer.Visible)));

%-------------------------------------------------
subplot(2,2,4);
hold on;
colors = lines(9);
hPP = gobjects(1,9);
for k = 1:9
    hPP(k) = plot(QB(:,k), PB(:,k), 'Color', colors(k,:));
end
hold off;
title('phase portrait blades');
xlabel('QB (u/m)');
ylabel('PB (Impulse / Geschwindigkeit)');
leg2 = legend(hPP, ...
    {'u_{11}','u_{12}','u_{13}','u_{21}','u_{22}','u_{23}','u_{31}','u_{32}','u_{33}'}, ...
    'Location','best');
grid on;
set(gcf,'Color','white');

% Interaktive Legende für Phase Portrait
set(leg2, 'ItemHitFcn', @(src, event) ...
    set(event.Peer, 'Visible', toggle(event.Peer.Visible)));

end

%-------------------------------------------------
% Sichtbarkeit umschalten
function newState = toggle(currentState)
    if strcmp(currentState, 'on')
        newState = 'off';
    else
        newState = 'on';
    end
end
