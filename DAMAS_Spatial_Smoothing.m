clear; close all; clc;
%==========================================================================
% DAMAS + Spatial Smoothing (SS) — SINAIS REAIS (WPE) do projeto IME
%
% Ideia: o array tem DUAS triades circulares iguais, deslocadas em altura:
%   sub-array 1 = mics [1 2 3] (z ~ 105 mm)
%   sub-array 2 = mics [4 5 6] (z ~ 265 mm)
% Trata cada triade como sub-array de 3 mics, calcula a covariancia 3x3 de
% cada uma e faz a MEDIA (forward spatial smoothing). Isso descorrelaciona
% (parcialmente) as fontes coerentes e devolve posto a Rx.
%
% Custo: a abertura efetiva cai de 6 -> 3 mics (resolucao menor).
% Aplica sobre o sinal WPE (mic*_wpe.wav).
%==========================================================================

c = 343; f = 700; M = 6;

%--------------------------------------------------------------
% Coordenadas dos mics (m)
%--------------------------------------------------------------
p1 = (1/100)*[237.0;   0.0;    105.89].';
p2 = (1/100)*[-118.60; 205.40; 105.89].';
p3 = (1/100)*[-118.60; -205.40;105.89].';
p4 = (1/100)*[-175.00; 0.0;    265.89].';
p5 = (1/100)*[87.50;  -151.60; 265.89].';
p6 = (1/100)*[87.50;   151.60; 265.89].';
P  = [p1; p2; p3; p4; p5; p6];

% sub-arrays (triades)
sub{1} = [1 2 3];
sub{2} = [4 5 6];
Msub = 3;                 % mics por sub-array

%--------------------------------------------------------------
% Le os 6 canais dereverberados (WPE)
%--------------------------------------------------------------
pasta = 'gravacoes/';
arqs = {'mic1_wpe.wav','mic2_wpe.wav','mic3_wpe.wav', ...
        'mic4_wpe.wav','mic5_wpe.wav','mic6_wpe.wav'};
[x1, fs] = audioread([pasta arqs{1}]);
tam = length(x1);
X = zeros(M, tam); X(1,:) = x1(:).';
for m = 2:M
    [xm, fsm] = audioread([pasta arqs{m}]);
    if fsm ~= fs, error('fs diferente entre canais.'); end
    L = min(tam, length(xm));
    X(m,1:L) = xm(1:L).';
end

%--------------------------------------------------------------
% Filtro BPF 700 Hz (zero delay) + analitico (Hilbert)
%--------------------------------------------------------------
fo = f; ORDEM = 1500; deltaf = 100;
freqs = [0 fo-deltaf fo-deltaf/2 fo+deltaf/2 fo+deltaf fs/2]/(fs/2);
ganhos = [0 0 1 1 0 0];
h = firpm(ORDEM, freqs, ganhos, [1 1 1]);
for m = 1:M
    xm = filtfilt(h,1,X(m,:).');
    X(m,:) = hilbert(xm).';
end

%==========================================================================
% SPATIAL SMOOTHING: media das covariancias 3x3 das duas triades
%==========================================================================
Rss = zeros(Msub, Msub);
for s = 1:numel(sub)
    Xs = X(sub{s}, :);              % 3 x tam
    Rs = (Xs*Xs')/tam;             % covariancia 3x3 da triade
    Rss = Rss + Rs;
end
Rss = Rss / numel(sub);            % media -> Rx suavizada (3x3)

%==========================================================================
% GRID
%==========================================================================
deltagraus = 1;
theta_vec = 60:deltagraus:120;
phi_vec   = 40:deltagraus:140;
Ntheta = numel(theta_vec); Nphi = numel(phi_vec);
N = Ntheta*Nphi;

%==========================================================================
% DAS com 3 mics (usa SO a 1a triade como referencia geometrica)
% steering de 3 mics; guarda p/ montar A
%==========================================================================
Pref = P(sub{1}, :);   % posicoes da triade de referencia (3x3)
Asteer = zeros(Msub, N);
b_img  = zeros(N,1);
col = 0;
for it = 1:Ntheta
    th = theta_vec(it);
    for ip = 1:Nphi
        ph = phi_vec(ip);
        col = col+1;
        vk = [sind(th)*cosd(ph); sind(th)*sind(ph); cosd(th)];
        a = exp(-1i*2*pi*f * Pref * vk / c);   % 3x1
        Asteer(:,col) = a;
        w = a/Msub;
        b_img(col) = real(w' * Rss * w);
    end
end
Imagem_DAS = reshape(b_img, Nphi, Ntheta).';

%==========================================================================
% Matriz A (PSF) com 3 mics  e  DAMAS
%==========================================================================
G = Asteer' * Asteer;
A = (1/Msub^2) * real(G .* conj(G));

NLOOPS = 500;
p = zeros(N,1); Adiag = diag(A);
disp('Rodando DAMAS (spatial smoothing)...'); tic
for loop = 1:NLOOPS
    for i = 1:N
        interf = A(i,:)*p - A(i,i)*p(i);
        p(i) = max(0, (b_img(i) - interf)/Adiag(i));
    end
end
fprintf('DAMAS+SS pronto em %d s.\n', round(toc));
Imagem_DAMAS = reshape(p, Nphi, Ntheta).';

%==========================================================================
% PLOTS
%==========================================================================
cA_phi=116.56; cA_theta=77.60; cB_phi=76.50; cB_theta=83.22;

figure
imagesc(phi_vec, theta_vec, Imagem_DAS/max(Imagem_DAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi (azimute)'); ylabel('\theta (zênite)')
title('DAS + Spatial Smoothing — sinais WPE')
hold on
plot(cA_phi,cA_theta,'w+','MarkerSize',15,'LineWidth',2)
plot(cB_phi,cB_theta,'wx','MarkerSize',15,'LineWidth',2)
legend('Caixa A','Caixa B')

figure
imagesc(phi_vec, theta_vec, Imagem_DAMAS/max(Imagem_DAMAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi (azimute)'); ylabel('\theta (zênite)')
title('DAMAS + Spatial Smoothing — sinais WPE')
hold on
plot(cA_phi,cA_theta,'w+','MarkerSize',15,'LineWidth',2)
plot(cB_phi,cB_theta,'wx','MarkerSize',15,'LineWidth',2)
legend('Caixa A','Caixa B')

%==========================================================================
% PICOS
%==========================================================================
[~,idx]=max(Imagem_DAS(:)); [tD,pD]=ind2sub(size(Imagem_DAS),idx);
fprintf('\n--- DAS+SS  pico:  theta=%.0f, phi=%.0f ---\n', theta_vec(tD), phi_vec(pD));
[~,idx]=max(Imagem_DAMAS(:)); [tM,pM]=ind2sub(size(Imagem_DAMAS),idx);
fprintf('--- DAMAS+SS pico: theta=%.0f, phi=%.0f ---\n', theta_vec(tM), phi_vec(pM));
fprintf('--- Caixa A: theta=%.1f, phi=%.1f ---\n', cA_theta, cA_phi);
fprintf('--- Caixa B: theta=%.1f, phi=%.1f ---\n', cB_theta, cB_phi);