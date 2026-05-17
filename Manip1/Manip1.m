clc; clear all; close all;

% load('manip1_y1_20kHz.mat');
load('manip1_y2_24kHz.mat');

freq_echantillonnage = 250000;
nb_echantillon = 400000;
% nb_echantillon=freq_echantillonnage*10;
t = (0:nb_echantillon-1)/freq_echantillonnage;

% figure;
% plot(t, y(:,2));
% hold on;
% plot(t, y(:,6), 'r');
% xlabel('Temps (s)');
% ylabel('Amplitude (v)');

% Retards mesurés à la main
T1 = 1.17e-3;
T2 = 1.27e-3;
T3 = 1.24e-3;
T4 = 1.37e-3;
T5 = 1.29e-3;
Delta_t = [T1; T2; T3; T4; T5];

% Positions des capteurs et emetteur en cm :
E  = [225, 100, 85];
R1 = [375, 100, 81];
R2 = [375, 200, 55];
R3 = [50,  50,  88];
R4 = [50, 120,  64];
R5 = [50, 200,  53];
R = [R1; R2; R3; R4; R5];

%% Calculer les Vi pour chaque capteur avec di/Ti
d1 = norm(E - R1) / 100;
d2 = norm(E - R2) / 100;
d3 = norm(E - R3) / 100;
d4 = norm(E - R4) / 100;
d5 = norm(E - R5) / 100;
d = zeros(5,1);
for i = 1:5
    d(i) = norm(E - R(i,:)) / 100; 
end

V1 = d1 / T1;
V2 = d2 / T2;
V3 = d3 / T3;
V4 = d4 / T4;
V5 = d5 / T5;

fprintf('Vitesse estimée V1 = %.2f m/s\n', V1);
fprintf('Vitesse estimée V2 = %.2f m/s\n', V2);
fprintf('Vitesse estimée V3 = %.2f m/s\n', V3);
fprintf('Vitesse estimée V4 = %.2f m/s\n', V4);
fprintf('Vitesse estimée V5 = %.2f m/s\n', V5);

V_moyenne = mean([V1 V2 V3 V4 V5]);
fprintf('Vitesse moyenne = %.2f m/s\n', V_moyenne);

%% Calculer les Vi pour chaque capteur avec argmin
V_est = (d.' * d) / (d.' * Delta_t);
% V_est = (Delta_t.' * d) / (Delta_t.' * Delta_t);
fprintf('Vitesse estimée = %.2f m/s\n', V_est);

%% Calculer les vitesses avec la corrélation des signaux Emetteur et Ri
signal_E = y(:,6);
signal_E = signal_E - mean(signal_E);
V = zeros(5,1);
Delta_t = zeros(5,1);

for i = 1:5
    signal_R = y(:,i);
    signal_R = signal_R - mean(signal_R);

    [c, lags] = xcorr(signal_R, signal_E, 'coeff');
    [~, idx] = max(abs(c));

    Delta_t(i) = abs(lags(idx)) / freq_echantillonnage;
    V(i) = d(i) / Delta_t(i);

    figure;
    plot(lags/freq_echantillonnage, c);
    xlabel('Temps (s)');
    ylabel('Corrélation');
    grid on;

    fprintf('R%d : Δt = %.3e s | V = %.2f m/s\n', i, Delta_t(i), V(i));
end

V_est = (d.' * d) / (d.' * Delta_t);
fprintf('Vitesse estimée (LS) = %.2f m/s\n', V_est);

d = [d1; d3; d4; d5];
D = [Delta_t(1); Delta_t(3:end)];
V_est = (d.' * d) / (d.' * D);
fprintf('Vitesse estimée (sans capteur 2) = %.2f m/s\n', V_est);


%%
% c = 1500;          % m/s
% R = R/100;         % m
% R = R.';           % 3 x N
% E = E(:);          % 3 x 1
% E = E/100;         % m
% 
% N = size(R,2);
% tu = zeros(N,1);
% 
% d = vecnorm(E - R, 2, 1);   % distances (m)
% tu = d.' / c;              % Nx1
% 
% disp('Temps de propagation théoriques (ms) :');
% disp(tu*1000);


%% Algorithme sphérique 
c = V_est;           % vitesse estimée
R = [R1; R3; R4; R5];
R = R/100; % en m
R = R.';             % 3 x N
E = E/100; % en m

N = size(R,2);

t_meas = [Delta_t(1); Delta_t(3:5)];
% t_meas = Delta_t;   % vecteur Nx1 des temps mesurés
M = mean(R,2);       % initialisation (centre capteurs)
 
epsilon = 1e-6;
maxIter = 20;
lambda = 1e-3;

for k = 1:maxIter

    d = vecnorm(M - R, 2, 1);     % distances Ri -> M
    r = d - c * t_meas.';         % résidus (1xN)

    J = zeros(N,3);
    for i = 1:N
        J(i,:) = (M - R(:,i)).'/d(i);
    end

    delta = -(J.' * J) \ (J.' * r.');
    % delta = -(J.'*J + lambda*eye(3)) \ (J.'*r.');

    M = M + delta;

    if norm(delta) < epsilon
        fprintf('Convergence atteinte à l''itération %d\n', k);
        break;
    end
end

fprintf('Position estimée de la source : [%.2f %.2f %.2f] m\n', M);
fprintf('Position reel de la source : [%.2f %.2f %.2f] m\n', E);
 
figure; hold on; grid on; axis equal;

scatter3(R(1,:), R(2,:), R(3,:), 100, 'b', 'filled');   % Capteurs
scatter3(E(1), E(2), E(3), 120, 'g', 'filled');         % Position vraie
scatter3(M(1), M(2), M(3), 120, 'r*');                  % Position estimée

scatter3(0, 0, 0, 150, '+k');                  % Origine (surface)

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Profondeur Z (m)');
title('Localisation sphérique');

legend({'Capteurs','Position vraie','Position estimée','Surface (0,0,0)'});

set(gca, 'ZDir', 'reverse');  
view(3);
