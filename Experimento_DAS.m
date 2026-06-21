clear; close all; clc;
%----------------
% Inicializações:
%----------------
c = 343; % velocidade do som em m/s
f = 700; % frequência do tom
M = 6;   % número de microfones
%---------------------------------
% Coordenadas dos mics (em metro):
%---------------------------------
p1 = (1/100)*[237.0; 0.0; 105.89].';
p2 = (1/100)*[-118.60; 205.40; 105.89].';
p3 = (1/100)*[-118.60; -205.40; 105.89].';
p4 = (1/100)*[-175.00; 0.0; 265.89].';
p5 = (1/100)*[87.50; -151.60; 265.89].';
p6 = (1/100)*[87.50; 151.60; 265.89].';
P  = [p1; p2; p3; p4; p5; p6]; % Matriz de posições (6x3)
%------------------------------
% Sinais reais (gravações):
%------------------------------
[x1, fs]  = audioread('gravacoes/mic1.wav');
[x2, fs2] = audioread('gravacoes/mic2.wav');
[x3, fs3] = audioread('gravacoes/mic3.wav');
[x4, fs4] = audioread('gravacoes/mic4.wav');
[x5, fs5] = audioread('gravacoes/mic5.wav');
[x6, fs6] = audioread('gravacoes/mic6.wav');

if ~isequal(fs, fs2, fs3, fs4, fs5, fs6)
    error('As taxas de amostragem das gravações não são iguais.');
end

% Truncar todos para o mesmo tamanho
tam = min([length(x1), length(x2), length(x3), length(x4), length(x5), length(x6)]);
x1 = x1(1:tam); x2 = x2(1:tam); x3 = x3(1:tam);
x4 = x4(1:tam); x5 = x5(1:tam); x6 = x6(1:tam);
%-----------------------------------------------------------------------
% Filtrando os sinais por um BPF centrado em 700 Hz (zero delay)
%-----------------------------------------------------------------------
fo = f;
ORDEM = 1500;
deltaf = 100;
freqs = [0 fo-deltaf fo-deltaf/2 fo+deltaf/2 fo+deltaf fs/2]/(fs/2);
ganhos = [0 0 1 1 0 0];
h = firpm(ORDEM, freqs, ganhos, [1 1 1]);

x1 = filtfilt(h,1,x1); x2 = filtfilt(h,1,x2);
x3 = filtfilt(h,1,x3); x4 = filtfilt(h,1,x4);
x5 = filtfilt(h,1,x5); x6 = filtfilt(h,1,x6);
%-----------------------------------------------------------------------
% Tornando os sinais analíticos e montando matrix snapshots
%-----------------------------------------------------------------------
x1 = hilbert(x1); x2 = hilbert(x2);
x3 = hilbert(x3); x4 = hilbert(x4);
x5 = hilbert(x5); x6 = hilbert(x6);
snapshots = [x1.'; x2.'; x3.'; x4.'; x5.'; x6.']; % 6 x tam
%-------------------------------------
% Estimando matriz de covariância Rx:
%-------------------------------------
Rx = snapshots * snapshots' / tam;
%-------------------------------------
% DAS (Delay-and-Sum):
%-------------------------------------
% O DAS é o beamformer mais simples:
% Aplica atrasos para alinhar os sinais de cada direção
% e soma — w = exp(-j*2*pi*f*atraso) / M
disp('Be patient and wait!')
tic

deltagraus = 1;
dimensaoy = length(0:deltagraus:180);
dimensaox = length(0:deltagraus:360);
Imagem = zeros(dimensaoy, dimensaox);

contaphi = 0;
for phi = 0:deltagraus:360
    contaphi = contaphi + 1;
    contatheta = 0;
    for theta = 0:deltagraus:180
        contatheta = contatheta + 1;

        % Vetor de direção para (theta, phi)
        vectork = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];

        % Vetor de pesos do DAS — atraso + normalização por M
        w = exp(-1i*2*pi*f * P * vectork / c) / M;

        % Potência na direção (theta, phi)
        Imagem(contatheta, contaphi) = real(w' * Rx * w);

    end
end

fprintf('You were patient enough to wait for %d seconds.\n', round(toc));

%-------------------------------------
% Plotando resultados
%-------------------------------------
phi_vec   = 0:deltagraus:360;
theta_vec = 0:deltagraus:180;

% Heatmap 2D
figure
imagesc(phi_vec, theta_vec, 10*log10(Imagem/max(Imagem(:))))
colorbar
colormap jet
xlabel('\phi (azimute em graus)')
ylabel('\theta (zênite em graus)')
title('DAS — Pseudoespectro (dB)')
hold on
plot(116.56, 77.60, 'w+', 'MarkerSize', 15, 'LineWidth', 2) % Caixa A
plot(76.50,  83.22,  'wx', 'MarkerSize', 15, 'LineWidth', 2) % Caixa B
legend('Caixa A esperada', 'Caixa B esperada')

% Heatmap 3D
figure
[PHI, THETA] = meshgrid(phi_vec, theta_vec);
surf(PHI, THETA, Imagem, 'EdgeColor', 'none')
axis vis3d; view(45,30); rotate3d on;
xlabel('\phi'); ylabel('\theta')
title('DAS — Heatmap 3D')

% Encontrar pico principal
[maxVal, idx] = max(Imagem(:));
[thetaIdx, phiIdx] = ind2sub(size(Imagem), idx);
fprintf('\n--- Pico encontrado pelo DAS ---\n')
fprintf('Pico principal: phi=%.1f graus, theta=%.1f graus\n', (phiIdx-1)*deltagraus, (thetaIdx-1)*deltagraus);
fprintf('Caixa A esperada: phi=116.56, theta=11.86\n');
fprintf('Caixa B esperada: phi=76.50, theta=6.27\n');