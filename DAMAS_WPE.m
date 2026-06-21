clear; close all; clc;
%==========================================================================
% DAMAS + FB (+ DL opcional) — SINAIS REAIS do projeto IME
% Pipeline DAMAS validado no simulado. Aqui tratamos a Rx ANTES de gerar o b:
%   FB = Forward-Backward Averaging -> quebra a correlacao das fontes
%   DL = Diagonal Loading (opcional) -> regulariza/estabiliza a Rx
% Objetivo: corrigir o deslocamento dos picos causado pelas fontes coerentes.
%==========================================================================

%----------------
% Inicializações:
%----------------
c = 343; % velocidade do som (m/s)
f = 700; % frequência do tom (Hz)
M = 6;   % número de microfones

%---------------------------------
% Coordenadas dos mics (em metro):
%---------------------------------
p1 = (1/100)*[237.0;   0.0;    105.89].';
p2 = (1/100)*[-118.60; 205.40; 105.89].';
p3 = (1/100)*[-118.60; -205.40;105.89].';
p4 = (1/100)*[-175.00; 0.0;    265.89].';
p5 = (1/100)*[87.50;  -151.60; 265.89].';
p6 = (1/100)*[87.50;   151.60; 265.89].';
P  = [p1; p2; p3; p4; p5; p6];  % Matriz de posições (6x3)

%==========================================================================
% (1) SINAIS REAIS (gravações)
%==========================================================================
[x1, fs]  = audioread('gravacoes/mic1.wav');
[x2, fs2] = audioread('gravacoes/mic2.wav');
[x3, fs3] = audioread('gravacoes/mic3.wav');
[x4, fs4] = audioread('gravacoes/mic4.wav');
[x5, fs5] = audioread('gravacoes/mic5.wav');
[x6, fs6] = audioread('gravacoes/mic6.wav');
% [x1, fs]  = audioread('gravacoes/mic1_wpe.wav');
% [x2, fs2] = audioread('gravacoes/mic2_wpe.wav');
% [x3, fs3] = audioread('gravacoes/mic3_wpe.wav');
% [x4, fs4] = audioread('gravacoes/mic4_wpe.wav');
% [x5, fs5] = audioread('gravacoes/mic5_wpe.wav');
% [x6, fs6] = audioread('gravacoes/mic6_wpe.wav');

if ~isequal(fs, fs2, fs3, fs4, fs5, fs6)
    error('As taxas de amostragem das gravações não são iguais.');
end

% Truncar todos para o mesmo tamanho
tam = min([length(x1), length(x2), length(x3), length(x4), length(x5), length(x6)]);
x1 = x1(1:tam); x2 = x2(1:tam); x3 = x3(1:tam);
x4 = x4(1:tam); x5 = x5(1:tam); x6 = x6(1:tam);

%-----------------------------------------------------------------------
% Filtrando por BPF centrado em 700 Hz (zero delay)
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
% Sinais analíticos + matriz de snapshots
%-----------------------------------------------------------------------
x1 = hilbert(x1); x2 = hilbert(x2);
x3 = hilbert(x3); x4 = hilbert(x4);
x5 = hilbert(x5); x6 = hilbert(x6);
snapshots = [x1.'; x2.'; x3.'; x4.'; x5.'; x6.']; % 6 x tam
K = tam;

%-------------------------------------
% Estimando matriz de covariância Rx:
%-------------------------------------
Rx = snapshots * snapshots' / K;

%==========================================================================
% (1b) TRATAMENTO DA Rx — Forward-Backward Averaging (+ Diagonal Loading)
%--------------------------------------------------------------------------
% FB: quebra a simetria persimétrica imposta pelas fontes correlacionadas.
%     Rx_fb = (Rx + J*conj(Rx)*J)/2,  J = matriz de troca (anti-identidade).
% DL: regulariza a Rx somando epsilon*I (epsilon = fração da potência média).
%     Mantido como FLAG para testar o efeito do FB isoladamente primeiro.
%==========================================================================
usar_FB = true;     % <-- Forward-Backward Averaging
usar_DL = false;    % <-- Diagonal Loading (ligue DEPOIS de avaliar o FB)
gamma_DL = 0.05;    % fração da potência média (traço/M) usada no loading

if usar_FB
    J = fliplr(eye(M));               % matriz de troca (exchange matrix)
    Rx = (Rx + J*conj(Rx)*J) / 2;     % Forward-Backward Averaging
    disp('FB aplicado na Rx.');
end

if usar_DL
    epsilon = gamma_DL * real(trace(Rx))/M;   % nível de loading
    Rx = Rx + epsilon*eye(M);                  % Diagonal Loading
    fprintf('DL aplicado na Rx (epsilon=%.3g).\n', epsilon);
end

%==========================================================================
% (2) GRID DE DIREÇÕES  -- largo o bastante p/ cobrir AS DUAS caixas
%     Caixa A: phi~116, theta~78  |  Caixa B: phi~76, theta~83
%==========================================================================
deltagraus = 1;                  % <-- 1 grau (sem gargalo). Depois: on-the-fly p/ 0.2
theta_vec = 60:deltagraus:120;   % zênite
phi_vec   = 40:deltagraus:140;   % azimute
Ntheta = numel(theta_vec);
Nphi   = numel(phi_vec);
N      = Ntheta*Nphi;
fprintf('Grid: %d direcoes (theta=%d, phi=%d). A sera %dx%d.\n', N, Ntheta, Nphi, N, N);

%==========================================================================
% (3) DAS  -- guarda steering vectors p/ reaproveitar no DAMAS
%==========================================================================
disp('Calculando DAS...'); tic
Asteer = zeros(M, N);
b_img  = zeros(N, 1);
col = 0;
for it = 1:Ntheta
    theta = theta_vec(it);
    for ip = 1:Nphi
        phi = phi_vec(ip);
        col = col + 1;
        vectork = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];
        a = exp(-1i*2*pi*f * P * vectork / c);
        Asteer(:,col) = a;
        w = a / M;
        b_img(col) = real(w' * Rx * w);
    end
end
fprintf('DAS pronto em %d s.\n', round(toc));
Imagem_DAS = reshape(b_img, Nphi, Ntheta).';

%==========================================================================
% (4) MATRIZ A (PSF)  -- Eq. (14)
%==========================================================================
disp('Montando matriz A (PSF)...'); tic
G = Asteer' * Asteer;
A = (1/M^2) * real(G .* conj(G));
fprintf('A montada (%dx%d) em %d s.\n', N, N, round(toc));

%==========================================================================
% (5) DAMAS  -- Algoritmo 1 / Eq. (16)
%==========================================================================
NLOOPS = 500;
p = zeros(N,1);
Adiag = diag(A);
disp('Rodando DAMAS...'); tic
for loop = 1:NLOOPS
    for i = 1:N
        interf = A(i,:)*p - A(i,i)*p(i);
        p(i) = max(0, (b_img(i) - interf)/Adiag(i));
    end
end
fprintf('DAMAS pronto (%d iteracoes) em %d s.\n', NLOOPS, round(toc));
Imagem_DAMAS = reshape(p, Nphi, Ntheta).';

%==========================================================================
% (6) PLOTS
%==========================================================================
% posições esperadas das caixas
cA_phi = 116.56; cA_theta = 77.60;
cB_phi = 76.50;  cB_theta = 83.22;

% string de status p/ títulos
status = 'DAMAS';
if usar_FB, status = [status ' + FB']; end
if usar_DL, status = [status ' + DL']; end

% ---- DAS 2D ----
figure
imagesc(phi_vec, theta_vec, Imagem_DAS/max(Imagem_DAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi (azimute)'); ylabel('\theta (zênite)')
title(['DAS — sinais reais (' status(7:end) ' na Rx)'])
hold on
plot(cA_phi, cA_theta, 'w+', 'MarkerSize',15,'LineWidth',2)
plot(cB_phi, cB_theta, 'wx', 'MarkerSize',15,'LineWidth',2)
legend('Caixa A','Caixa B')

% ---- DAMAS 2D ----
figure
imagesc(phi_vec, theta_vec, Imagem_DAMAS/max(Imagem_DAMAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi (azimute)'); ylabel('\theta (zênite)')
title([status ' — sinais reais'])
hold on
plot(cA_phi, cA_theta, 'w+', 'MarkerSize',15,'LineWidth',2)
plot(cB_phi, cB_theta, 'wx', 'MarkerSize',15,'LineWidth',2)
legend('Caixa A','Caixa B')

% ---- DAMAS 3D ----
figure
[PHI,THETA] = meshgrid(phi_vec, theta_vec);
surf(PHI, THETA, Imagem_DAMAS/max(Imagem_DAMAS(:)), 'EdgeColor','none')
view(45,30); xlabel('\phi'); ylabel('\theta'); zlabel('Normalizado')
title([status ' — 3D (sinais reais)'])

%==========================================================================
% (7) PICOS
%==========================================================================
[~, idx] = max(Imagem_DAS(:));
[tD,pD] = ind2sub(size(Imagem_DAS), idx);
fprintf('\n--- DAS  pico:   theta=%.0f, phi=%.0f ---\n', theta_vec(tD), phi_vec(pD));
[~, idx] = max(Imagem_DAMAS(:));
[tM,pM] = ind2sub(size(Imagem_DAMAS), idx);
fprintf('--- DAMAS pico:  theta=%.0f, phi=%.0f ---\n', theta_vec(tM), phi_vec(pM));
fprintf('--- Caixa A esperada: theta=%.1f, phi=%.1f ---\n', cA_theta, cA_phi);
fprintf('--- Caixa B esperada: theta=%.1f, phi=%.1f ---\n', cB_theta, cB_phi);