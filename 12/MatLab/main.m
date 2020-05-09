% Inggeo Uebung 12
% 13.Mai 2020
% Ziqing Yu 3218051
clc
clear all
close all
%% Import Data
load data.mat 

%% Aufgabe a
xi_1 = data (1:20,4) - data(1:20,3); % Höhenanomalie 1 - 20
sigma_HN = 0.001; % Standardabweichung Normalhöhen
sigma_e = 0.005;  % Standardabweichung Ellipslid Höhe
sigma_xi = sqrt(sigma_HN^2 + sigma_e^2);  % Fehlerfortpflanzung
%
figure
plot(xi_1)

%% Aufgabe b
A_1 = [ones(20,1), data(1:20,1), data(1:20,2), data(1:20,1).* data(1:20,2), data(1:20,1).^2, data(1:20,2).^2];  % Matrix bauen
a_list = (A_1' * A_1) \ A_1' * xi_1;  % Ausgleichen

%% Aufgabe c

%% Aufgabe d
A_2 = [ones(10,1), data(21:30,1), data(21:30,2), data(21:30,1).* data(21:30,2), data(21:30,1).^2, data(21:30,2).^2]; 
xi_2 = A_2 * a_list;    % Höhenanomalie 21 - 30
NH_under =  data(21:30,4) - xi_2; % Normalhöhen 21 - 30
data(21:30,3) = NH_under;

%% Aufgabe e
