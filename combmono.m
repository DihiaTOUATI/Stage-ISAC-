clc; clear; close all; rng(0);
c = 3e8;        % vitesse de la lumière (m/s)
N = 64;         % nombre de sous-porteuses OFDM
delta_f = 15e3; % espacement fréquentiel entre sous-porteuses
fs = N * delta_f; % fréquence d'échantillonnage garantissant l'orthogonalité OFDM
Nsym = 128;      % nombre de symboles OFDM (dimension temporelle pour Doppler)
f0 = 10e9;       % fréquence porteuse radar
lambda = c / f0; % longueur d'onde associée
B = 20e6;        % bande passante du chirp
T = 1e-4;        % durée du chirp
t = (0:1/fs:T-1/fs)'; % vecteur temps discret
Nt = length(t);       % nombre d'échantillons du chirp
k = B / T;            % pente du chirp (variation de fréquence)

chirp = exp(1j*2*pi*(0.5*k*t.^2)); % génération du chirp LFM (fréquence croissante)

bits = randi([0 1], 2*N, Nsym);  % génération aléatoire de bits (QPSK : 2 bits/symbole)
X_all = zeros(N, Nsym);          % matrice contenant les symboles OFDM (fréquence)

for m = 1:Nsym      % boucle sur les symboles OFDM
    data = bi2de(reshape(bits(:,m), [], 2)); % conversion bits → symboles (0 à 3)
    X = exp(1j*(pi/2*data));        % mapping QPSK : % 0 → 1   % 1 → j  % 2 → -1   % 3 → -j
    X(1) = 1;                       % insertion d'un pilote connu
    X_all(:,m) = X;                 % on stocke chaque symboles OFDM dans une colone m
end

x_ofdm = ifft(X_all);          % passage au domaine temporel (signal OFDM)

x_all = zeros(N, Nsym);        % matrice qui va contenir le signal combiné OFDM + chirp

for m = 1:Nsym                 % combinaison OFDM + chirp pour chaque symbole
    x_all(:,m) = x_ofdm(:,m) .* chirp(1:N); % On combine OFDM et chirp en multipliant les deux signaux.
                       % on met l'information OFDM sur one onde radar 
end

d = 1000;    % distance réelle de la cible
v = 30;      % vitesse réelle de la cible
tau = 2*d/c; % délai aller-retour du signal
f_D = 2*v / lambda; % fréquence Doppler due au mouvement
delay_samples = round(tau * fs); % conversion du délai en nombre d'échantillons

y_all = zeros(size(x_all)); %  matrice pour stocker le signal recu 

for m = 1:Nsym      % boucle canal pour chaque symbole
    h = zeros(1, delay_samples+1); % creer un  vecteur et on place 1 a la position du retard 
    h(delay_samples+1) = 1;  % crée un retard pur
    y = conv(x_all(:,m), h); % convolution → applique le retard au signal 
    y = y(1:N);              % tronque pour garder la taille originale
    y = y .* exp(1j*2*pi*f_D*m*(N/fs)); % applique le Doppler (variation de phase lente)
    y_all(:,m) = y; % stockage du signal reçu
end

% ===================== DISTANCE =====================
[corr, lags] = xcorr(y_all(:,1), x_all(:,1));  % corrélation pour estimer le retard
[~, idx] = max(abs(corr));  % détection du pic de corrélation
tau_est = lags(idx)/fs;     % conversion du lag en temps
d_est = c * tau_est / 2;    % conversion temps → distance

% ===================== Doppler  =====================
doppler_fft = fftshift(fft(y_all, [], 2),2); % FFT sur les symboles (analyse Doppler)
f_doppler = (-Nsym/2:Nsym/2-1)/(Nsym*(N/fs)); % axe des fréquences Doppler
power = sum(abs(doppler_fft).^2,1); % puissance Doppler (somme sur sous-porteuses)
[~, idx_d] = max(power); % détection du pic Doppler
f_D_est = f_doppler(idx_d); % fréquence Doppler estimée
v_est = (f_D_est * lambda) / 2; % conversion Doppler → vitesse

fprintf('Distance vraie = %.2f m\n', d); 
fprintf('Distance estimée = %.2f m\n\n', d_est);
fprintf('Vitesse vraie = %.2f m/s\n', v); 
fprintf('Doppler estimé = %.2f Hz\n', f_D_est); 
fprintf('Vitesse estimée = %.2f m/s\n', v_est); 
figure; 

subplot(2,1,1); 
plot(lags/fs, abs(corr));
title('Corrélation (Distance)'); 
grid on;

subplot(2,1,2); 
plot(f_doppler, power); 
title('Spectre Doppler'); 
xlabel('Fréquence (Hz)');
grid on; 