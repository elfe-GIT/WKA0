function [t,y] = solveIVP(sys,caseID, initID)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% set-upy numerical paramters
sys.i = caseID;
if sys.Oga(caseID) == 0
    sys.tEnd = 5;
else
    % one full rotation
    sys.tEnd = 2*pi()/sys.Oga(caseID);
end
tspan = [0,sys.tEnd];

% initial values from sys structured for ODE45/ODE15s
y0          = zeros(28,1);
y0(initID,1)= 1;
%
h = waitbar(0,'solving ODE - please wait...');
[t,y] = ode15s(@(t,y)wkadydt(t,y,sys),tspan,y0);
close(h);

end