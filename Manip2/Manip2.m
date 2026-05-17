clc; clear; close all;

% load('manip2_y2_stastion.mat');   % stastionnaire
% nb_echantillon = 750000;

load('manip2_y1_nonstastion.mat');   % non-stastionnaire
nb_echantillon = 1250000;

fs = 250000;          % Hz
c  = 1314.46;            % m/s
t = (0:nb_echantillon-1)/fs;

% figure;
% plot(t,y(:,1));
% xlabel('Temps (s)');
% ylabel('Amplitude (v)');

% Positions des capteurs (cm)
R = [375 100 81;
     375 200 55;
      50  50 88;
      50 120 64;
      50 200 53];

R = R / 100;          % m
R = R.';              % 3 x N
N = size(R,2);

E  = [275, 40, 20];
E = E /100;

Delta_t = zeros(N-1,1);

signal_ref = y(:,1);
signal_ref = signal_ref - mean(signal_ref);

for i = 2:N
    signal_i = y(:,i);
    signal_i = signal_i - mean(signal_i);

    maxLag = round(0.0001 * fs);
    % maxLag = round(0.0005 * fs);
    [corr, lags] = xcorr(signal_i, signal_ref, maxLag, 'coeff');

    [~, idx] = max(abs(corr));
    Delta_t(i-1) = lags(idx) / fs;

    % fprintf('Δt(1-%d) = %.3e s\n', i, Delta_t(i-1));

    % figure;
    % plot(lags/fs, corr);
    % xlabel('Temps (s)');
    % ylabel('Corrélation');
    % grid on;
end

% %% TDOA théoriques (pour validation)
% 
% Delta_t_th = zeros(N-1,1);
% 
% d1 = norm(E - R(:,1));
% 
% for i = 2:N
%     di = norm(E - R(:,i));
%     Delta_t_th(i-1) = (di - d1) / c;
% end
% 
% disp('--- Comparaison TDOA ---')
% for i = 2:N
%     fprintf('Capteur %d : Δt_mes = %.3e s | Δt_th = %.3e s\n', i, Delta_t(i-1), Delta_t_th(i-1));
% end

% Initialisation (centre des capteurs)
M = mean(R,2);
epsilon = 1e-6;
maxIter = 30;

for k = 1:maxIter

    d1 = norm(M - R(:,1));
    r  = zeros(N-1,1);
    J  = zeros(N-1,3);

    for i = 2:N
        di = norm(M - R(:,i));

        % Résidu hyperbolique
        r(i-1) = (di - d1) - c * Delta_t(i-1);

        % Jacobien
        ei = (M - R(:,i)) / di;
        e1 = (M - R(:,1)) / d1;

        J(i-1,:) = (ei - e1).';
    end

    % Correction
    delta = -(J.'*J) \ (J.'*r);
    % delta = - (J.'*J + 1e-6*eye(3)) \ (J.'*r);

    M = M + delta;

    if norm(delta) < epsilon
        fprintf('Convergence atteinte à l’’itération %d\n', k);
        break;
    end
end

fprintf('\nPosition estimée du bateau : [%.2f %.2f %.2f] m\n', M);
fprintf('Position reel du bateau : [%.2f %.2f %.2f] m\n', E);

figure; hold on; grid on; axis equal;
scatter3(R(1,:), R(2,:), R(3,:), 100, 'b', 'filled');
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
