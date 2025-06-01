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
sys.i = 2;
% set-upy numerical paramters
sys.tEnd = 2*pi()/sys.Oga(sys.i);
sys.tEnd = 1.
tspan = [0,sys.tEnd];
% initial values from sys structured for ODE45
y0     = zeros(28,1);
y0(5,1)=1;
%
h = waitbar(0,'solving ODE - please wait...');
[t,y] = ode15s(@(t,y)wkadydt(t,y,sys),tspan,y0);
close(h);
%% postprocess
% define plotting range
range = [0.,1];
postprocess(sys,t,y,range);

%% EOF

