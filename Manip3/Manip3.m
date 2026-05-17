clc; clear; close all;

load('manip3_GH013152.mat');   % y : Ns x 5

fs = 250000;          
c  = 1314.46;            

[Ns, ~] = size(y);

t = (0:Ns-1)/fs;

figure;
plot(t,y);
xlabel('Temps (s)');
ylabel('Amplitude (v)');

%% =========================
% Marquer visuellement entrée et sortie
%% =========================
t_entree = 4.11;   % temps d'entrée du bateau (en s)
t_sortie = 8;   % temps de sortie du bateau (en s)

%% =========================
% Découper y pour ne garder que la zone d’intérêt
%% =========================
idx_entree = round(t_entree*fs) + 1;
idx_sortie = round(t_sortie*fs);

y = y(idx_entree:idx_sortie, :);
t = t(idx_entree:idx_sortie);

[Ns, N] = size(y);

%% =========================
% Positions des capteurs (m)
%% =========================
R = [375 100 81;
     375 200 55;
      50  50 88;
      50 120 64;
      50 200 53];
R = R / 100;      
R = R.';          % 3 x 5

%% =========================
% Paramètres temporels
%% =========================
window_length = 0.02;              % 20 ms
step = 0.01;                       % recouvrement 50 %
L = round(window_length * fs);
step_samples = round(step * fs);
nb_windows = floor((Ns - L)/step_samples);

E_est = zeros(nb_windows,3);
time_vec = zeros(nb_windows,1);

%% =========================
% Boucle temporelle
%% =========================
for k = 1:nb_windows
    
    idx_start = (k-1)*step_samples + 1;
    idx_end   = idx_start + L - 1;
    segment = y(idx_start:idx_end,:);
    signal_ref = segment(:,1);
    signal_ref = signal_ref - mean(signal_ref);
    Delta_t = zeros(N-1,1);
    
    % === Estimation TDOA ===
    for i = 2:N
        signal_i = segment(:,i);
        signal_i = signal_i - mean(signal_i);
        maxLag = round(0.003 * fs);   % ~750 samples
        [corr,lags] = xcorr(signal_i,signal_ref,maxLag,'coeff');
        [~,idx] = max(abs(corr));
        Delta_t(i-1) = lags(idx)/fs;
    end
    
    %% =========================
    % Estimation position instantanée
    %% =========================
    
    M = mean(R,2);      % initialisation
    epsilon = 1e-6;
    maxIter = 20;
    
    for iter = 1:maxIter
        
        d1 = norm(M - R(:,1));
        r  = zeros(N-1,1);
        J  = zeros(N-1,3);
        
        for i = 2:N
            di = norm(M - R(:,i));
            
            r(i-1) = (di - d1) - c*Delta_t(i-1);
            
            ei = (M - R(:,i))/di;
            e1 = (M - R(:,1))/d1;
            
            J(i-1,:) = (ei - e1).';
        end
        
        delta = -(J.'*J + 1e-6*eye(3)) \ (J.'*r);
        M = M + delta;
        
        if norm(delta) < epsilon
            fprintf('Convergence atteinte à l’itération %d\n', iter);
            break
        end
    end
    
    E_est(k,:) = M.';
    time_vec(k) = mean(t(idx_start:idx_end));
    
end

Px = polyfit(time_vec, E_est(:,1), 1);
Py = polyfit(time_vec, E_est(:,2), 1);
Pz = polyfit(time_vec, E_est(:,3), 1);

x_fit = polyval(Px, time_vec);
y_fit = polyval(Py, time_vec);
z_fit = polyval(Pz, time_vec);

% idx = (x_fit >= 0.3) & (x_fit <= 4);
% x_fit = x_fit(idx);
% y_fit = y_fit(idx);

%% =========================
% Load video
%% =========================
videoFile = 'Project.mp4';
v = VideoReader(videoFile);
fps = v.FrameRate;
T_total = v.Duration;

%% =========================
% Known sensor positions (meters)
%    (X Y Z) -> we only use X,Y
%% =========================
world_pts = R(1:2, :).'; % use only X,Y (planar assumption)

%% =========================
% Read FIRST frame
%% =========================
frame1 = readFrame(v);

figure;
imshow(frame1);
hold on;
title('Cliquez les 5 capteurs dans l''ordre');

x_pix = zeros(5,1);
y_pix = zeros(5,1);

for i = 1:5
    [x_pix(i), y_pix(i)] = ginput(1);
    
    % Affichage du point rouge
    plot(x_pix(i), y_pix(i), 'ro', 'MarkerSize',10, 'LineWidth',2);
    
    % Afficher le numéro du capteur
    text(x_pix(i)+5, y_pix(i), sprintf('%d', i), ...
        'Color','red','FontSize',30,'FontWeight','bold');
end

image_pts = [x_pix y_pix];

%% =========================
% Compute homography (image -> world)
%% =========================
tform = fitgeotrans(image_pts, world_pts, 'projective');

%% =========================
% Click boat INITIAL position
%% =========================
figure;
imshow(frame1);
title('Cliquez sur la position initiale du bateau');

[x0_pix, y0_pix] = ginput(1);
[x0, y0] = transformPointsForward(tform, x0_pix, y0_pix);

%% =========================
% Read LAST frame
%% =========================
v.CurrentTime = v.Duration - 1/fps;
frameLast = readFrame(v);

figure;
imshow(frameLast);
title('Cliquez sur la position finale du bateau');

[x1_pix, y1_pix] = ginput(1);
[x1, y1] = transformPointsForward(tform, x1_pix, y1_pix);

%% =========================
% Generate straight-line trajectory
%% =========================
N = 200;  % number of trajectory points

t_real = linspace(0, T_total, N);

x_real = x0 + (x1 - x0) * (t_real / T_total);
y_real = y0 + (y1 - y0) * (t_real / T_total);

%% =========================
% Affichage trajectoire
%% =========================

figure; hold on; grid on; axis equal;

scatter3(R(1,:),R(2,:),R(3,:),100,'b','filled');
% plot(E_est(:,1), E_est(:,2), 'r*');
plot(x_fit, y_fit, 'k','LineWidth', 2);
plot(x_real, y_real, 'r','LineWidth',2);

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Trajectoire estimée');
legend('Capteurs','Trajectoire estimée ajustée', 'Trajectoire réel');
xlim([-1, 6])
set(gca,'ZDir','reverse');
view(3);

