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
% Usar 25 - 30 segundos
%inicio = 25 * fs;
%fim = 30 * fs;
%x1 = x1(inicio:fim); x2 = x2(inicio:fim); x3 = x3(inicio:fim);
%x4 = x4(inicio:fim); x5 = x5(inicio:fim); x6 = x6(inicio:fim);
%-----------------------------------------------------------------------
% Filtrando os sinais por um BPF centrado em 700 Hz (zero delay)
%-----------------------------------------------------------------------
fo = f;
ORDEM = 1500;
deltaf = 100;
freqs = [0 fo-deltaf fo-deltaf/2 fo+deltaf/2 fo+deltaf fs/2]/(fs/2);
ganhos = [0 0 1 1 0 0];
h = firpm(ORDEM, freqs, ganhos, [1 1 1]);

x1 = filtfilt(h,1,x1);
x2 = filtfilt(h,1,x2);
x3 = filtfilt(h,1,x3);
x4 = filtfilt(h,1,x4);
x5 = filtfilt(h,1,x5);
x6 = filtfilt(h,1,x6);
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
% MUSIC:
%-------------------------------------
% Número de fontes (suas 2 caixas de som)
nFontes = 2;

% Decomposição em autovalores
[E, D] = eig(Rx);
autovalores = diag(D);

% Ordenar autovalores do maior para o menor
[~, idx] = sort(autovalores, 'descend');
E = E(:, idx); % Reordena autovetores

% Subespaço do ruído = autovetores associados aos (M - nFontes) menores autovalores
En = E(:, nFontes+1:end); % 6x4 — os 4 autovetores de ruído

disp('Be patient and wait!')
tic

% Grid azimute (phi de 0 a 360) e zênite (theta de 0 a 180)
deltagraus = 1; % Use 1 grau para ser mais rápido; mude para 0.1 para mais precisão
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

        % Vetor de steering (direção de busca)
        vaux = exp(-1i*2*pi*f * P * vectork / c); % 6x1

        % Pseudoespectro MUSIC:
        % Quanto menor a projeção de vaux no subespaço do ruído,
        % maior o pico — indica presença de fonte nessa direção
        ruido = vaux' * (En * En') * vaux;
        Imagem(contatheta, contaphi) = 1 / abs(ruido);

    end
end

clc
fprintf('You were patient enough to wait for %d seconds.\n', round(toc));

%-------------------------------------
% Plotando resultados
%-------------------------------------
phi_vec   = 0:deltagraus:360;
theta_vec = 0:deltagraus:180;

% Heatmap 2D (mais fácil de analisar)
figure
imagesc(phi_vec, theta_vec, 10*log10(Imagem/max(Imagem(:))))
colorbar
colormap jet
xlabel('\phi (azimute em graus)')
ylabel('\theta (zênite em graus)')
title('MUSIC — Pseudoespectro (dB)')
% Marcando os DOAs esperados das caixas
hold on
plot(116.56, 77.60, 'w+', 'MarkerSize', 15, 'LineWidth', 2) % Caixa A (zenite correto)
plot(76.02,  83.22, 'wx', 'MarkerSize', 15, 'LineWidth', 2) % Caixa B (zenite correto)
legend('Caixa A esperada', 'Caixa B esperada')

% Heatmap 3D
figure
[PHI, THETA] = meshgrid(phi_vec, theta_vec);
surf(PHI, THETA, Imagem, 'EdgeColor', 'none')
axis vis3d; view(45,30); rotate3d on;
xlabel('\phi'); ylabel('\theta')
title('MUSIC — Heatmap 3D')

% Encontrar e exibir os 2 picos principais
Imagem_temp = Imagem;
fprintf('\n--- Picos encontrados pelo MUSIC ---\n')
for k = 1:2
    [~, idx] = max(Imagem_temp(:));
    [tIdx, pIdx] = ind2sub(size(Imagem_temp), idx);
    fprintf('Pico %d: phi=%.1f graus, theta=%.1f graus\n', k, (pIdx-1)*deltagraus, (tIdx-1)*deltagraus);
    % Suprimir região ao redor do pico para achar o próximo
    Imagem_temp(max(1,tIdx-10):min(end,tIdx+10), max(1,pIdx-10):min(end,pIdx+10)) = 0;
end