clc
close all
clear all

TS16GR4 = importfile16("E:\Studium\5,6-Ingenieurgeodaesie\Uebung\IngGeo-6-Semester-Uebung\14\matlab\TS16_GR4_mod.txt", [1, Inf]);
TS30GR4 = importfile30("E:\Studium\5,6-Ingenieurgeodaesie\Uebung\IngGeo-6-Semester-Uebung\14\matlab\TS30_GR4_mod.txt", [1, Inf]);

R = 6378137;

z1 = zeros(10,1);
z1(1:5) = TS16GR4(11:15,5) / 200 * pi;
z1(6:10) = (400 - TS16GR4(16:20,5)) / 200 * pi;

z2 = zeros(10,1);
z2(1:5) = TS30GR4(11:15,5) / 200 * pi;
z2(6:10) = (400 - TS30GR4(16:20,5)) / 200 * pi;

sr1 = TS16GR4(11:20,6);
sr2 = TS30GR4(11:20,6);

sh1 = TS16GR4(11:20,7);
sh2 = TS30GR4(11:20,7);
sh = mean([sh1;sh2]);
dH = (sr1 .* cos(z1) - sr2 .* cos(z2))/2;
dH_mean = mean(dH);

% Fehlerfortpflanzung
F = zeros(10,40);
diag_sigma_zs = zeros(1,40);
for i = 1:10
    F(i,4 * i - 3) = cos(z1(i)) / 2;
    F(i,4 * i - 2) = sr1(i) / 2 * (-sin(z1(i)));
    F(i,4 * i - 1) = -cos(z2(i)) / 2;
    F(i,4 * i - 0) = sr2(i) / 2 * (sin(z2(i)));
    
    diag_sigma_zs(4 * i -3) = 1e-3 + 1.5e-6 * z1(i);
    diag_sigma_zs(4 * i -2) = 0.15e-3 / 200 * pi;
    diag_sigma_zs(4 * i -1) = 1e-3 + 1e-6 * z2(i);
    diag_sigma_zs(4 * i -0) = 0.3e-3 / 200 * pi;
end
Sigma_zs_quad = diag(diag_sigma_zs.^2);
Sigma_h_quad = F * Sigma_zs_quad * F';
sigma_h = sqrt(diag(Sigma_h_quad));

% eine Vereinfachung
z1q = zeros(10,1);
z1q(1:5) = TS16GR4(1:5,5) / 200 * pi;
z1q(6:10) = (400 - TS16GR4(6:10,5)) / 200 * pi;

z2q = zeros(10,1);
z2q(1:5) = TS30GR4(1:5,5) / 200 * pi;
z2q(6:10) = (400 - TS30GR4(6:10,5)) / 200 * pi;
k = 1 + (200 - z1q/pi*200 - z2q/pi*200) / ((200 / pi) .* (sh / R));
k_mean = mean(k);

% fuer Vergleich
dh_ue10 = 481.2389-474.2540;