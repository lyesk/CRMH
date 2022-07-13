%% Convert .csv reaction rates to .mat collision frequencies
% Source of .csv files: http://www.eirene.de/juniplib/eiram/
% AMJUEL H.4 database
% Data was downloaded by specifying electron density ne=1e12 cm-3 for all reactions
% Roughly 10000 points between 0.1 and 100 eV electron temperature

ne=1e12;
%% Input .csv data:
%First column is electron temperature in eV
%Second column is collision rate in cm^3/s

inputFolder='/home/lyes/Documents/PhD/CRM_H/Eirene_data/csv/';
filename='H.4_2.2.12';
T=readtable([inputFolder filename '.csv']);

Te=table2array(T(:,1));
coll_rate=table2array(T(:,2));

%% Output .mat data:
%First column is electron temperature in eV
%Second column is collision rate in cm^3/s

outputArray=[Te coll_rate];

outputFolder='/home/lyes/Documents/PhD/CRM_H/Eirene_data/mat/';
outputFilename=string([outputFolder filename '.mat']);

save(outputFilename,'outputArray');