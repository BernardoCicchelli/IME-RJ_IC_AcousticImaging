clear; close all; clc;
%==========================================================================
% DAMAS — Reprodução do exemplo simulado do Tech Report (Apolinario, 2026)
% Reproduz a Fig. 2 (DAS) e a Fig. 3 (DAMAS) com sinal SINTÉTICO.
% Geometria: 6 mics reais do projeto IME. Grid grosso (1 grau) p/ validar.
%==========================================================================

%----------------
% Inicializações:
%----------------
c  = 343;    % velocidade do som (m/s)
f  = 700;    % frequência do tom (Hz)
M  = 6;      % número de microfones
fs = 44100;  % taxa de amostragem (Hz)

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
% (1) GERAÇÃO DO SINAL SIMULADO  -- substitui o bloco de audioread
%--------------------------------------------------------------------------
% Reproduz a Sec. III do paper: WGN passa-baixa modulado por portadora de
% 700 Hz, vindo de UMA fonte em (theta_SOI, phi_SOI), com ruído branco
% espacialmente descorrelacionado para SNR = 10 dB.
% Para voltar aos sinais REAIS: comente este bloco e descomente o audioread.
%==========================================================================
theta_SOI = 95;    % zênite da fonte (graus)
phi_SOI   = 252;   % azimute da fonte (graus)
SNR_dB    = 10;    % relação sinal-ruído (dB)
T         = 15;    % duração (s)  -> paper usa 15 s (661.500 amostras)
K         = round(T*fs);  % número de snapshots

wo = 2*pi*f/fs;    % frequência digital da portadora

% --- sinal-fonte: WGN passa-baixa (banda-base ~ deltaf) ---
deltaf = 100;                            % largura de banda do envelope (Hz)
s = randn(1, K);                         % WGN
[bb, ab] = butter(4, (deltaf/2)/(fs/2)); % passa-baixa
s = filter(bb, ab, s);                   % envelope band-limited
s = s / sqrt(mean(abs(s).^2));           % normaliza potência

% --- vetor de direção e steering da fonte ---
vectork_SOI = [sind(theta_SOI)*cosd(phi_SOI); ...
               sind(theta_SOI)*sind(phi_SOI); ...
               cosd(theta_SOI)];
% atraso geométrico (em amostras) de cada mic
tau = (P * vectork_SOI) / c;             % atraso em segundos (6x1)

% --- monta sinais analíticos de cada mic: s(k) modulado + fase do atraso ---
k_idx = 0:K-1;
carrier = exp(1i*wo*k_idx);              % portadora analítica
snapshots = zeros(M, K);
for m = 1:M
    phase_delay = exp(-1i*2*pi*f*tau(m));   % fase do atraso (Eq.4)
    snapshots(m,:) = s .* carrier * phase_delay;
end

% --- adiciona ruído branco espacial p/ SNR alvo ---
Psig = mean(abs(snapshots(:)).^2);
Pn   = Psig / (10^(SNR_dB/10));
noise = sqrt(Pn/2)*(randn(M,K) + 1i*randn(M,K));
snapshots = snapshots + noise;

%==========================================================================
% (1-REAL) BLOCO ALTERNATIVO PARA SINAIS REAIS (deixar comentado por ora)
%==========================================================================
% [x1, fs]  = audioread('gravacoes/mic1.wav');
% ... (seu bloco original de audioread + filtro FIR + hilbert) ...
% snapshots = [x1.'; x2.'; ... ; x6.'];
%==========================================================================

%-------------------------------------
% Estimando matriz de covariância Rx:
%-------------------------------------
Rx = snapshots * snapshots' / K;

%==========================================================================
% (2) GRID DE DIREÇÕES  (exemplo do paper, mas com passo grosso p/ validar)
%==========================================================================
deltagraus = 1;                 % <-- comece com 1 grau; depois 0.2
theta_vec = 70:deltagraus:110;  % zênite
phi_vec   = 210:deltagraus:270; % azimute
Ntheta = numel(theta_vec);
Nphi   = numel(phi_vec);
N      = Ntheta*Nphi;           % total de direções no grid

%==========================================================================
% (3) DAS  -- reproduz a Fig. 2
% Também guardamos TODOS os steering vectors p/ reaproveitar no DAMAS.
%==========================================================================
disp('Calculando DAS...'); tic

Asteer = zeros(M, N);   % matriz com todos os steering vectors (M x N)
b_img  = zeros(N, 1);   % imagem DAS empilhada (vetor coluna) -> "b"
col = 0;
for it = 1:Ntheta
    theta = theta_vec(it);
    for ip = 1:Nphi
        phi = phi_vec(ip);
        col = col + 1;
        vectork = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];
        a = exp(-1i*2*pi*f * P * vectork / c);   % steering (NÃO normalizado)
        Asteer(:,col) = a;
        w = a / M;                               % peso DAS (Eq.10)
        b_img(col) = real(w' * Rx * w);          % potência (Eq.11)
    end
end
fprintf('DAS pronto em %d s.\n', round(toc));

% reorganiza b em matriz (theta x phi) para visualizar
Imagem_DAS = reshape(b_img, Nphi, Ntheta).';

%==========================================================================
% (4) MATRIZ A  (Point-Spread Function)  -- Eq. (14)
% A = (1/M^2) * (Asteer^H Asteer) .* conj(Asteer^H Asteer)
%==========================================================================
disp('Montando matriz A (PSF)...'); tic
G = Asteer' * Asteer;            % N x N (produtos internos entre direções)
A = (1/M^2) * (G .* conj(G));    % Hadamard -> |.|^2  (real, não-negativa)
A = real(A);
fprintf('A montada (%dx%d) em %d s.\n', N, N, round(toc));

%==========================================================================
% (5) DAMAS  -- Algoritmo 1 / Eq. (16)
% Gauss-Seidel com projeção em não-negatividade.
%==========================================================================
NLOOPS = 500;
p = zeros(N,1);                  % imagem deconvoluída (inicia em zero)
Adiag = diag(A);

disp('Rodando DAMAS...'); tic
for loop = 1:NLOOPS
    for i = 1:N
        % soma dos outros termos (j != i) usando p atual (Gauss-Seidel)
        interf = A(i,:)*p - A(i,i)*p(i);
        p(i) = max(0, (b_img(i) - interf)/Adiag(i));
    end
end
fprintf('DAMAS pronto (%d iteracoes) em %d s.\n', NLOOPS, round(toc));

Imagem_DAMAS = reshape(p, Nphi, Ntheta).';

%==========================================================================
% (6) PLOTS  -- Figs. 2 e 3
%==========================================================================
% ---- DAS 2D ----
figure
imagesc(phi_vec, theta_vec, Imagem_DAS/max(Imagem_DAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi (azimute)'); ylabel('\theta (zênite)')
title('DAS — imagem acústica (Fig. 2)')
hold on; plot(phi_SOI, theta_SOI, 'w+', 'MarkerSize',15,'LineWidth',2)
legend('DoA correto')

% ---- DAS 3D ----
figure
[PHI,THETA] = meshgrid(phi_vec, theta_vec);
surf(PHI, THETA, Imagem_DAS/max(Imagem_DAS(:)), 'EdgeColor','none')
view(45,30); xlabel('\phi'); ylabel('\theta'); zlabel('Normalizado')
title('DAS — 3D')

% ---- DAMAS 2D ----
figure
imagesc(phi_vec, theta_vec, Imagem_DAMAS/max(Imagem_DAMAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi (azimute)'); ylabel('\theta (zênite)')
title('DAMAS — imagem acústica (Fig. 3)')
hold on; plot(phi_SOI, theta_SOI, 'w+', 'MarkerSize',15,'LineWidth',2)
legend('DoA correto')

% ---- DAMAS 3D ----
figure
surf(PHI, THETA, Imagem_DAMAS/max(Imagem_DAMAS(:)), 'EdgeColor','none')
view(45,30); xlabel('\phi'); ylabel('\theta'); zlabel('Normalizado')
title('DAMAS — 3D')

%==========================================================================
% (7) PICOS
%==========================================================================
[~, idx] = max(Imagem_DAS(:));
[tD,pD] = ind2sub(size(Imagem_DAS), idx);
fprintf('\n--- DAS  pico:   theta=%.0f, phi=%.0f ---\n', theta_vec(tD), phi_vec(pD));
[~, idx] = max(Imagem_DAMAS(:));
[tM,pM] = ind2sub(size(Imagem_DAMAS), idx);
fprintf('--- DAMAS pico:  theta=%.0f, phi=%.0f ---\n', theta_vec(tM), phi_vec(pM));
fprintf('--- Esperado:    theta=%d, phi=%d ---\n', theta_SOI, phi_SOI);