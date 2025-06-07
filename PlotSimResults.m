% PlotSimResults.m
% Lädt gespeicherte Simulationsergebnisse und führt den Postprozess aus

% Datei laden aus dem Unterordner "Results"
load(fullfile('Results', 'SimResu_Parm1_Case13_4.mat'), 't', 'y', 'sys');

% Zeitbereich für Plot definieren
range = [0, 0];  % [0,0] = gesamte Zeit

% Postprocessing
postprocessTWO(sys, t, y, range);
