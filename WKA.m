%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Matlab-Script                                                     %
% Uses Version: R2021a                                              %
% Modul Maschinendynamik                                            %
% Löst: Das Anfangswertproblem für die WKA                          %
% Autor: A. Baumgart                                                %
% Last Update: 2025-05-27                                           %
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cd('C:\Users\wzo448\OneDrive - HAW-HH\2025-SS\Maschinendynamik\Matlab\WKA0'); 
addpath("./Classes","./Functions");

% start .....
parameterFile='parameter.xlsx';
sheetName='Systemparameter';
%% preprocess
sys = preprocess(parameterFile,sheetName);

%% solve IVP
sys.i = 5;
% set-upy numerical paramters
if sys.Oga(sys.i) == 0
    sys.tEnd = 5;
else
    sys.tEnd = 2*pi()/sys.Oga(sys.i);
end

tspan = [0,sys.tEnd];
% initial values from sys structured for ODE45
y0      = zeros(28,1);
y0( 1,1)= 1;
%
h = waitbar(0,'solving ODE - please wait...');
[t,y] = ode15s(@(t,y)wkadydt(t,y,sys),tspan,y0);
close(h);
%% postprocess
% define plotting range
range = [0.,sys.tEnd];
postprocess(sys,t,y,range);

%% now: characteristic multiplyers

% characteristic multiplyers
C = zeros(28,length(sys.Oga));
for i = 1:length(sys.Oga)
    X = sprintf('finished %3.0f %%', i/length(sys.Oga)*100);
    % select Omega
    sys.i = i;
    % set-upy numerical paramters
    if sys.Oga(sys.i) ~= 0
        sys.tEnd = 2*pi()/sys.Oga(sys.i);
        tspan = [0,sys.tEnd];
        % monodromy matrix
        MDM = zeros(28,28);
        for j=1:28
            % initial values from sys structured for ODE45
            y0      = zeros(28,1);
            y0( j,1)= 1;
            %
            h = waitbar(0,'solving ODE - please wait...');
            [t,y] = ode15s(@(t,y)wkadydt(t,y,sys),tspan,y0);
            close(h);
            MDM(:,j) = transpose(y(end,:));
        end
        C(:,i) = eig(MDM);
        disp(X);
    end
end


%% EOF
