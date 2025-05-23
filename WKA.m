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

%% postprocess
% define plotting range
postprocess(sys,t,y,range);

%% EOF

