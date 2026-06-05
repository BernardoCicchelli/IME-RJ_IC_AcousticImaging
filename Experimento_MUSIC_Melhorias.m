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
% MELHORIA 0c — Ruído branco filtrado (sinal descorrelacionado)
% Em vez de usar os .wav diretamente, adicionamos ruído branco filtrado
% para descorrelacionar as fontes e testar o sistema
%-----------------------------------------------------------------------
% Gera ruído branco com o mesmo tamanho dos sinais
rng(42); % semente fixa para reprodutibilidade
ruido1 = randn(tam, 1);
ruido2 = randn(tam, 1);
ruido3 = randn(tam, 1);
ruido4 = randn(tam, 1);
ruido5 = randn(tam, 1);
ruido6 = randn(tam, 1);

% Filtra o ruído na mesma banda de 700 Hz
ruido1 = filtfilt(h,1,ruido1); ruido2 = filtfilt(h,1,ruido2);
ruido3 = filtfilt(h,1,ruido3); ruido4 = filtfilt(h,1,ruido4);
ruido5 = filtfilt(h,1,ruido5); ruido6 = filtfilt(h,1,ruido6);

% Mistura os sinais reais com o ruído branco filtrado
% O ruído descorrelaciona as fontes mantendo a banda de interesse
alpha = 0.3; % fator de mistura — ajuste entre 0 (só real) e 1 (só ruído)
x1 = x1 + alpha * ruido1; x2 = x2 + alpha * ruido2;
x3 = x3 + alpha * ruido3; x4 = x4 + alpha * ruido4;
x5 = x5 + alpha * ruido5; x6 = x6 + alpha * ruido6;

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
%-----------------------------------------------------------------------
% MELHORIA 1b — Janela curta e média por segmentos
% Calcula Rx em janelas de 256ms e faz média temporal
% Atenua componentes tardios da reverberação
%-----------------------------------------------------------------------
%faz um for comecando em 10ms e indo ate 500ms, depois plotando com 2
%autovalores
janela_ms = 256; % tamanho da janela em ms
janela_amostras = round(janela_ms * fs / 1000); % converte para amostras
num_janelas = floor(tam / janela_amostras); % número de janelas

Rx = zeros(M, M); % inicializa Rx zerada

for j = 1:num_janelas
    % Pega o trecho de cada janela
    inicio_j = (j-1) * janela_amostras + 1;
    fim_j = j * janela_amostras;

    % Snapshots da janela atual
    snap_j = snapshots(:, inicio_j:fim_j);

    % Soma a contribuição dessa janela
    Rx = Rx + snap_j * snap_j' / janela_amostras;
end

% Média sobre todas as janelas
Rx = Rx / num_janelas;
% Forward-Backward Averaging
%J = fliplr(eye(M));
%Rx = (Rx + J * conj(Rx) * J) / 2;

% Diagonal Loading
epsilon = 0.01 * trace(Rx) / M;
Rx = Rx + epsilon * eye(M);
%-------------------------------------
% MUSIC:
%-------------------------------------
nFontes = 2;

[E, D] = eig(Rx);
autovalores = diag(D);
[autoval_sorted, idx] = sort(autovalores, 'descend');
E = E(:, idx);

fprintf('Autovalores ordenados:\n')
disp(real(autoval_sorted))

En = E(:, nFontes+1:end); % subespaço do ruído

disp('Be patient and wait!')
tic

% MELHORIA 0b — Grid limitado a theta < 120° (reduz falsos picos)
deltagraus = 1;
theta_max = 120; % antes era 180°
dimensaoy = length(0:deltagraus:theta_max);
dimensaox = length(0:deltagraus:360);
Imagem = zeros(dimensaoy, dimensaox);

contaphi = 0;
for phi = 0:deltagraus:360
    contaphi = contaphi + 1;
    contatheta = 0;
    for theta = 0:deltagraus:theta_max % MELHORIA 0b aqui
        contatheta = contatheta + 1;

        vectork = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];
        vaux = exp(-1i*2*pi*f * P * vectork / c);
        ruido = vaux' * (En * En') * vaux;
        Imagem(contatheta, contaphi) = 1 / abs(ruido);

    end
end

fprintf('You were patient enough to wait for %d seconds.\n', round(toc));

%-------------------------------------
% Plotando resultados
%-------------------------------------
phi_vec   = 0:deltagraus:360;
theta_vec = 0:deltagraus:theta_max;

% Heatmap 2D
figure
imagesc(phi_vec, theta_vec, 10*log10(Imagem/max(Imagem(:))))
colorbar
colormap jet
xlabel('\phi (azimute em graus)')
ylabel('\theta (zênite em graus)')
title('MUSIC + 0b + 0c — Pseudoespectro (dB)')
hold on
plot(116.56, 77.60, 'w+', 'MarkerSize', 15, 'LineWidth', 2) % Caixa A
plot(76.02,  83.22, 'wx', 'MarkerSize', 15, 'LineWidth', 2) % Caixa B
legend('Caixa A esperada', 'Caixa B esperada')

% Heatmap 3D
figure
[PHI, THETA] = meshgrid(phi_vec, theta_vec);
surf(PHI, THETA, Imagem, 'EdgeColor', 'none')
axis vis3d; view(45,30); rotate3d on;
xlabel('\phi'); ylabel('\theta')
title('MUSIC + 0b + 0c — Heatmap 3D')

% Picos principais
Imagem_temp = Imagem;
fprintf('\n--- Picos encontrados pelo MUSIC + 0b + 0c ---\n')
for k = 1:2
    [~, idx] = max(Imagem_temp(:));
    [tIdx, pIdx] = ind2sub(size(Imagem_temp), idx);
    fprintf('Pico %d: phi=%.1f graus, theta=%.1f graus\n', k, (pIdx-1)*deltagraus, (tIdx-1)*deltagraus);
    Imagem_temp(max(1,tIdx-10):min(end,tIdx+10), max(1,pIdx-10):min(end,pIdx+10)) = 0;
end