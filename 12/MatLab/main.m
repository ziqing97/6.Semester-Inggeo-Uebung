% Inggeo Uebung 12
% 13.Mai 2020
% Ziqing Yu 3218051
clc
clear all
close all
%% Import Data
load data.mat

%% Aufgabe 1
xi = data (1:20,4) - data(1:20,3); % Höhenanomalie 1 - 20
sigma_HN = 0.001; % Standardabweichung Normalhöhen
sigma_e = 0.005;  % Standardabweichung Ellipslid Höhe
sigma_xi = sqrt(sigma_HN^2 + sigma_e^2);  % Fehlerfortpflanzung
%
figure
plot(xi)

%% Aufgabe 2 
