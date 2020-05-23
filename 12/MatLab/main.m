% Inggeo Uebung 12
% 13.Mai 2020
% Ziqing Yu 3218051
clc
clear all
close all
%% Import Data
load data.mat 

%% Aufgabe a
zeta_1 = data(1:20,4) - data(1:20,3); % Höhenanomalie 1 - 20
sigma_HN = 0.001; % Standardabweichung Normalhöhen
sigma_e = 0.005;  % Standardabweichung Ellipslid Höhe
sigma_zeta = sqrt(sigma_HN^2 + sigma_e^2);  % Fehlerfortpflanzung

xq = min(data(1:20,2)):50:max(data(1:20,2));
yq = min(data(1:20,1)):50:max(data(1:20,1));
[xq,yq] = meshgrid(xq,yq);
vq = griddata(data(1:20,2),data(1:20,1),zeta_1,xq,yq);
figure,hold on 
scatter3(data(1:20,2),data(1:20,1),zeta_1)
mesh(xq,yq,vq)
xlabel("x")
ylabel("y")
zlabel("Höhenanomalie")


%% Aufgabe b
A_1 = [ones(20,1), data(1:20,1), data(1:20,2), data(1:20,1).* data(1:20,2), data(1:20,1).^2, data(1:20,2).^2];  % Matrix bauen
a_bar = (A_1' * A_1) \ A_1' * zeta_1;  % Ausgleichen

a_bar = inv(A_1' * A_1) * A_1' * zeta_1; % test



%% Aufgabe c
r = 20 - length(a_bar);
zeta_1_bar = A_1 * a_bar;
epsilon_bar = zeta_1 - zeta_1_bar;
sigma_zeta_bar = sqrt(epsilon_bar' * epsilon_bar) / r;
Sigma_a = sigma_zeta_bar * inv(A_1' * A_1);

% Sigma_a = sigma_zeta * inv(A_1' * A_1); % test


sigma_a = sqrt(diag(Sigma_a));
T = abs(a_bar - 0) ./ sigma_a;
Q = [2.447, 2.571, 2.776, 3.182, 4.303, 12.71];   % Quantil
idx = find(T < Q(1));


a_list = cell(6,1);
T_list = cell(6,1);

%%
i = 1;
id = zeros(6,1) * NaN;

check = zeros(6,1) * NaN;
check_list = 1:6;

while ~isempty(idx)
    a_list{i} = a_bar;
    T_list{i} = T;
    id(i) = find(T == min(T));
    check(i) = check_list(id(i));
    check_list(id(i)) = [];
    A_1(:,id(i)) = [];
    a_bar = (A_1' * A_1) \ A_1' * zeta_1;
    
    a_bar = inv(A_1' * A_1) * A_1' * zeta_1; % test
    
    r = 20 - length(a_bar);
    zeta_1_bar = A_1 * a_bar;
    epsilon_bar = zeta_1 - zeta_1_bar;
    sigma_zeta_bar = sqrt(epsilon_bar' * epsilon_bar) / r;
    Sigma_a = sigma_zeta_bar^2 * inv(A_1' * A_1);
    
%     Sigma_a = sigma_zeta * inv(A_1' * A_1); % test
    
    sigma_a = sqrt(diag(Sigma_a));
    T = abs(a_bar - 0) ./ sigma_a;
    idx = find(T < Q(i+1));
    i = i + 1;
end

check = sort(check);  % Welche Elemente sind gelöscht

%% Aufgabe d
A_2 = [ones(10,1), data(21:30,1), data(21:30,2), data(21:30,1).* data(21:30,2), data(21:30,1).^2, data(21:30,2).^2]; 
for i=1:6
    if isnan(id(i)) 
        break
    else
        A_2(:,id(i)) = [];
    end
end
zeta_2 = A_2 * a_bar;    % Höhenanomalie 21 - 30
NH_under =  data(21:30,4) - zeta_2; % Normalhöhen 21 - 30
data(21:30,3) = NH_under;


%% Aufgabe e
F = [eye(10),A_2];
[~,l] = size(F);
Sigma_big = zeros(l,l);
Sigma_big(1:10,1:10) = 0.005^2 * eye(10);
Sigma_big(11:l,11:l) = Sigma_a;

Sigma_nh = F * Sigma_big * F';

sigma_nh = sqrt(diag(Sigma_nh));

