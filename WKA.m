%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Matlab-Script                                                     %
% Uses Version: R2021a                                              %
% Modul Maschinendynamik                                            %
% Löst: Das Anfangswertproblem für die WKA                          %
% Autor: A. Baumgart                                                %
% Last Update: 2025-05-27                                           %
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cd('C:\Users\JanAt\Documents\GitHub\WKA0'); 
addpath("./Classes","./Functions");
tic;
% start .....
parameterFile='parameter.xlsx';
sheetName='Systemparameter';
%% preprocess
sys = preprocess(parameterFile,sheetName);

%% solve IVP
caseID = 13; % select Omega-ID
initID = 3; % if of state variable with initial deflection 

cID = caseID; % for setting the filename at save results
iID = initID;
[t,y] = solveIVP(sys,caseID,initID);

%% postprocess
% define plotting range
range = [0.,t(end)];
postprocess(sys,t,y,range);

%% now: characteristic multiplyers

% characteristic multiplyers
C = zeros(28,length(sys.Oga));
for caseID = 1:length(sys.Oga)
    % set-upy numerical paramters
    if sys.Oga(caseID) ~= 0
        X = sprintf('finished %3.0f %%', caseID/length(sys.Oga)*100);
        % monodromy matrix
        MDM = zeros(28,28);
        for initID=1:28
            % initial values from sys structured for ODE45
            [t,y] = solveIVP(sys,caseID,initID);   
            MDM(:,initID) = transpose(y(end,:));
            fprintf('.');
        end
        C(:,caseID) = eig(MDM);
        disp(X);
    end
end

save(strcat("C-",string(datetime('now','Format','yyyy-MM-dd')),".mat"),"C");

% save results in folder "Results"
 outDir = 'Results';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
filename = fullfile(outDir, sprintf('SimResu_Parm1_Case_%d_%d.mat', cID, iID));
save(filename, 't', 'y', 'sys');

elapsedTime = toc/60
% and plot ....
hold off;
plot(sys.Oga,abs(C(:,:)),'*');
title('characteristic multiplyers')
xlabel('f/Hz ->');
ylabel('ρ/1 ->');

%% EOF
