clc; clear; close all; rng(0); 

c = 3e8;                  % Vitesse de la lumière (m/s)
N = 64;                   % nombre de sous-porteuses OFDM (taille FFT/IFFT)
delta_f = 15e3;           % espacement fréquentiel entre sous-porteuses (Hz)
fs = N * delta_f;         % fréquence d'échantillonnage du signal,  garantit l'orthogonalite 
Nsym = 32;                % nombre de symboles OFDM , 32 symboles OFDM donc on observe 32 instants differents 
f0 = 10e9;                % fréquence porteuse 
lambda = c / f0;          % longueur d'onde associée à la fréquence porteuse

% ===================== [BLOC 1] BIT STREAM =====================
bits = randi([0 1], 2*N, Nsym);      % génération aléatoire de bits (0 ou 1)
% bits = [0,1,1,0....]

% ===================== [BLOC 2] MODULATEUR =====================
X_all = zeros(N, Nsym);    % creer une matrice vide de dimension 64*32 pour stocker les symboles OFDM  

for m = 1:Nsym
    
    data = bi2de(reshape(bits(:,m), [], 2));    % on regroupe 2 bits pour chaque symbole 
                                                % convertit chaque paire en entier (0 à 3)

    X = exp(1j * (pi/2 * data));        % transformer les bits en signal complexe 
                                        % mapping QPSK : % 0 → 1   % 1 → j  % 2 → -1   % 3 → -j

    % ===================== [BLOC 4] PILOT =====================
    X(1) = 1;  % on force la 1 sous porteuse a etre connu pour estimer le canal 
                                        

    X_all(:,m) = X;                     % on stocke chaque symboles OFDM dans une colone m
end

% ===================== [BLOC 5] IFFT =====================
x_all = ifft(X_all);                   % transforme le signal du domaine fréquentiel vers le domaine temporel
                                       % crée le signal OFDM temporel

% ===================== [BLOC 7] CP =====================
d = 1000;                              % distance réelle de la cible (m)
v = 30;                                % vitesse de la cible (m/s)
tau = 2*d/c;                           % délai aller-retour du signal
delay_samples = round(tau * fs);       % On convertit le délai en secondes en nombre d'échantillons pour pouvoir appliquer le retard dans le signal discret
cp = max(ceil(0.1*N), delay_samples);  % longueur du préfixe cyclique , doit être au moins supérieur au retard pour éviter l'ISI
x_cp = [x_all(end-cp+1:end,:); x_all]; % ajout du préfixe cyclique : on copie la fin du signal et on la met au début

% ===================== [BLOC 9] CHANNEL =====================
f_D = 2*v / lambda;                    % fréquence Doppler due au mouvement de la cible
y_all = zeros(size(x_cp));             % matrice remplie 0 pour stocker les signaux reçus meme taille que x_cp
for m = 1:Nsym                         % création du canal 
    h = zeros(1, delay_samples+1);     % creer un  vecteur et on place 1 a la position du retard 
    h(delay_samples+1) = 1;            
    y = conv(x_cp(:,m), h);            % convolution → applique le retard au signal 
                                       
    y = y .* exp(1j*2*pi*f_D*m*(N/fs));  % ajoute un décalage de phase progressif pour modéliser l'effet Doppler
    y_all(:,m) = y(1:size(x_cp,1));    % on prends seulement une partie de y car la convolution agrandit le signal
end
% ===================== [BLOC 11] REMOVE CP =====================
rx = y_all(cp+1:cp+N, :);              % suppression du préfixe cyclique
                                       % récupération du signal utile

% ===================== [BLOC 13] FFT =====================
Y = fft(rx);                           % retour au domaine fréquentiel , récupération des sous-porteuses modifier par le canal

% ===================== [BLOC 15] EQUALIZATION =====================

H = Y ./ (X_all + 1e-10);              % estimation du canal , H = Y / X

% ===================== [BLOC 16] IFFT =====================
h_est = ifft(H, N);                    % transformation vers domaine temporel car h est en frequence et retard en temps 
% ===================== DISTANCE =====================
h_mean = mean(abs(h_est),2);           % moyenne sur les symboles OFDM , améliore la robustesse

[~, idx] = max(h_mean);                % trouve position du  pic (retard)
tau_est = (idx-1)/fs;                  % conversion indice → temps
d_est = c * tau_est / 2;               % conversion temps → distance

% ===================== DOPPLER =====================
doppler_fft = fftshift(fft(H, [], 2),2);       % FFT sur les colonne on analyse evolution entre symbole OFDM , permet d'extraire le Doppler
f_doppler = (-Nsym/2:Nsym/2-1)/(Nsym*(N/fs));  % creer l'axe des fréquences  Doppler
power = sum(abs(doppler_fft),1);               % puissance Doppler (somme sur sous-porteuses)
[~, idx_d] = max(power);                       % détection du pic Doppler
f_D_est = f_doppler(idx_d);                    % fréquence Doppler estimée
v_est = (f_D_est * lambda) / 2;                % conversion Doppler → vitesse

% ===================== AFFICHAGE =====================
fprintf('Distance vraie = %.2f m\n', d);
fprintf('Distance estimée = %.2f m\n\n', d_est);

fprintf('Vitesse vraie = %.2f m/s\n', v);
fprintf('Doppler estimé = %.2f Hz\n', f_D_est);
fprintf('Vitesse estimée = %.2f m/s\n', v_est);

% ===================== PLOTS =====================
figure;

subplot(2,1,1);
stem(h_mean,'filled');  
% affiche le profil de retard
title('Profil de retard OFDM');

subplot(2,1,2);
plot(f_doppler, power);  
% affiche le spectre Doppler
title('Spectre Doppler OFDM');
xlabel('Fréquence (Hz)');
grid on;