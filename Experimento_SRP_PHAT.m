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
% SRP-PHAT — Steered Response Power with Phase Transform
%-----------------------------------------------------------------------
% Diferença do DAS:
% Em vez de usar a potência do sinal diretamente,
% o SRP-PHAT normaliza o espectro cruzado pela sua amplitude
% mantendo APENAS a fase — isso torna o algoritmo robusto
% a reflexões e fontes correlacionadas
%
% Para cada par de microfones (l,k):
% GCC-PHAT(tau) = IFFT( X_l(w) * conj(X_k(w)) / |X_l(w) * conj(X_k(w))| )
%
% O SRP-PHAT soma as GCC-PHAT de todos os pares
% avaliadas no atraso esperado para cada direção
%-----------------------------------------------------------------------
% Todos os pares de microfones
pares = nchoosek(1:M, 2);
nPares = size(pares, 1);
fprintf('Número de pares de microfones: %d\n', nPares);
% Usa NFFT menor para melhor estabilidade
NFFT = 2^nextpow2(tam/100);
fprintf('NFFT usado: %d\n', NFFT);

% Divide o sinal em segmentos e faz média das GCC-PHAT
n_seg = floor(tam / NFFT);
GCC = zeros(NFFT, nPares);

for seg = 1:n_seg
    ini = (seg-1)*NFFT + 1;
    fim = seg*NFFT;

    % FFT de cada segmento
    X1s = fft(x1(ini:fim), NFFT); X2s = fft(x2(ini:fim), NFFT);
    X3s = fft(x3(ini:fim), NFFT); X4s = fft(x4(ini:fim), NFFT);
    X5s = fft(x5(ini:fim), NFFT); X6s = fft(x6(ini:fim), NFFT);
    Xs = [X1s, X2s, X3s, X4s, X5s, X6s];

    for p = 1:nPares
        l = pares(p,1); k = pares(p,2);
        Xlk = Xs(:,l) .* conj(Xs(:,k));
        Xlk_phat = Xlk ./ (abs(Xlk) + 1e-10);
        GCC(:,p) = GCC(:,p) + real(ifft(Xlk_phat, NFFT));
    end
end

% Média sobre os segmentos
GCC = GCC / n_seg;

disp('Be patient and wait!')
tic

% Grid limitado a theta < 120° (melhoria 0b)
deltagraus = 1;
theta_max = 120;
dimensaoy = length(0:deltagraus:theta_max);
dimensaox = length(0:deltagraus:360);
Imagem = zeros(dimensaoy, dimensaox);

contaphi = 0;
for phi = 0:deltagraus:360
    contaphi = contaphi + 1;
    contatheta = 0;
    for theta = 0:deltagraus:theta_max
        contatheta = contatheta + 1;

        % Vetor de direção
        vectork = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];

        % Para cada par, calcula o atraso esperado e busca na GCC-PHAT
        srp = 0;
        for p = 1:nPares
            l = pares(p,1);
            k = pares(p,2);

            % Atraso esperado entre microfones l e k para essa direção
            % tau = (distância do mic l à fonte - distância do mic k à fonte) / c
            tau = (P(l,:) - P(k,:)) * vectork / c; % em segundos

            % Converte atraso para índice na GCC-PHAT
            idx_tau = round(tau * fs) + 1; % +1 porque MATLAB começa em 1

            % Garante que o índice está dentro dos limites
            if idx_tau < 1
                idx_tau = idx_tau + NFFT;
            elseif idx_tau > NFFT
                idx_tau = idx_tau - NFFT;
            end

            % Soma a GCC-PHAT no atraso esperado
            srp = srp + GCC(idx_tau, p);
        end

        Imagem(contatheta, contaphi) = srp;
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
imagesc(phi_vec, theta_vec, Imagem)
colorbar
colormap jet
xlabel('\phi (azimute em graus)')
ylabel('\theta (zênite em graus)')
title('SRP-PHAT — Pseudoespectro')
hold on
plot(116.56, 77.60, 'w+', 'MarkerSize', 15, 'LineWidth', 2) % Caixa A (zênite correto)
plot(76.02,  83.22, 'wx', 'MarkerSize', 15, 'LineWidth', 2) % Caixa B (zênite correto)
legend('Caixa A esperada', 'Caixa B esperada')

% Heatmap 3D
figure
[PHI, THETA] = meshgrid(phi_vec, theta_vec);
surf(PHI, THETA, Imagem, 'EdgeColor', 'none')
axis vis3d; view(45,30); rotate3d on;
xlabel('\phi'); ylabel('\theta')
title('SRP-PHAT — Heatmap 3D')

% Picos principais
Imagem_temp = Imagem;
fprintf('\n--- Picos encontrados pelo SRP-PHAT ---\n')
for k = 1:2
    [~, idx] = max(Imagem_temp(:));
    [tIdx, pIdx] = ind2sub(size(Imagem_temp), idx);
    fprintf('Pico %d: phi=%.1f graus, theta=%.1f graus\n', k, (pIdx-1)*deltagraus, (tIdx-1)*deltagraus);
    Imagem_temp(max(1,tIdx-10):min(end,tIdx+10), max(1,pIdx-10):min(end,pIdx+10)) = 0;
end

fprintf('\nCaixa A esperada: phi=116.56, theta=77.60\n');
fprintf('Caixa B esperada: phi=76.02, theta=83.22\n');