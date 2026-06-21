clear; close all; clc;
%==========================================================================
% DAMAS on-the-fly (Delta=0.2 grau) — mapa largo, sinais WPE, 6 mics
% Ataca o GARGALO COMPUTACIONAL: a matriz A (N x N) seria ~190 GB em 0.2 grau.
% Solucao (do paper): NUNCA montar A inteira. Guardar so a Asteer (6 x N) e
% gerar cada LINHA de A on-the-fly, vetorizada, dentro do loop do DAMAS.
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

%--------------------------------------------------------------
% Le os 6 canais dereverberados (WPE) + BPF 700 Hz + Hilbert
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
    L = min(tam, length(xm)); X(m,1:L) = xm(1:L).';
end
fo=f; ORDEM=1500; deltaf=100;
freqs=[0 fo-deltaf fo-deltaf/2 fo+deltaf/2 fo+deltaf fs/2]/(fs/2);
ganhos=[0 0 1 1 0 0];
h=firpm(ORDEM,freqs,ganhos,[1 1 1]);
for m=1:M
    xm=filtfilt(h,1,X(m,:).');
    X(m,:)=hilbert(xm).';
end
Rx = (X*X')/tam;

%==========================================================================
% GRID FINO (0.2 grau) — mapa largo
%==========================================================================
deltagraus = 0.2;
theta_vec = 60:deltagraus:120;
phi_vec   = 40:deltagraus:140;
Ntheta=numel(theta_vec); Nphi=numel(phi_vec);
N = Ntheta*Nphi;
fprintf('Grid 0.2 grau: N=%d direcoes. (A inteira seria %.1f GB)\n', ...
        N, (N^2*8)/1e9);

%==========================================================================
% Pre-computa Asteer (6 x N) e b (DAS) — cabe na RAM (~MB)
%==========================================================================
disp('Pre-computando steering vectors e b (DAS)...'); tic
Asteer = zeros(M, N);
col=0;
for it=1:Ntheta
    th=theta_vec(it); st=sind(th); ct=cosd(th);
    for ip=1:Nphi
        ph=phi_vec(ip); col=col+1;
        vk=[st*cosd(ph); st*sind(ph); ct];
        Asteer(:,col)=exp(-1i*2*pi*f * P * vk / c);
    end
end
% b_i = a_i^H Rx a_i / M^2  (vetorizado)
RA = Rx*Asteer;                          % 6 x N
b_img = real(sum(conj(Asteer).*RA,1)).'/ (M^2);   % N x 1
fprintf('Pronto em %d s.\n', round(toc));

% diagonal de A: |a_i^H a_i|^2/M^2 = (M^2)/M^2 = 1 (steering nao normalizado)
% mas calculamos para garantir:
Adiag = real(sum(conj(Asteer).*Asteer,1)).'.^2 / (M^2);  % N x 1

%==========================================================================
% DAMAS on-the-fly — gera cada linha de A dentro do loop, vetorizada
%==========================================================================
NLOOPS = 2;      % menos iteracoes (paper: >100 pouco ajuda e custa caro)
p = zeros(N,1);
fprintf('Rodando DAMAS on-the-fly (%d iter, N=%d)... pode demorar.\n', NLOOPS, N);
tic
for loop=1:NLOOPS
    for i=1:N
        % linha i de A, on-the-fly:  |a_i^H * Asteer|^2 / M^2
        g = Asteer(:,i)' * Asteer;          % 1 x N (produtos internos)
        Ai = (abs(g).^2) / (M^2);           % 1 x N (linha i de A)
        interf = Ai*p - Ai(i)*p(i);
        p(i) = max(0, (b_img(i) - interf)/Adiag(i));
    end
    if mod(loop,10)==0
        fprintf('  iter %d/%d  (%.0f s decorridos)\n', loop, NLOOPS, toc);
    end
end
fprintf('DAMAS on-the-fly concluido em %d s.\n', round(toc));

Imagem_DAS   = reshape(b_img, Nphi, Ntheta).';
Imagem_DAMAS = reshape(p,     Nphi, Ntheta).';

%==========================================================================
% PLOTS
%==========================================================================
cA_phi=116.56; cA_theta=77.60; cB_phi=76.50; cB_theta=83.22;

figure
imagesc(phi_vec, theta_vec, Imagem_DAS/max(Imagem_DAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi'); ylabel('\theta'); title('DAS 0.2° — WPE')
hold on; plot(cA_phi,cA_theta,'w+',cB_phi,cB_theta,'wx','MarkerSize',15,'LineWidth',2)

figure
imagesc(phi_vec, theta_vec, Imagem_DAMAS/max(Imagem_DAMAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi'); ylabel('\theta'); title('DAMAS on-the-fly 0.2° — WPE')
hold on; plot(cA_phi,cA_theta,'w+',cB_phi,cB_theta,'wx','MarkerSize',15,'LineWidth',2)

[~,idx]=max(Imagem_DAMAS(:)); [tM,pM]=ind2sub(size(Imagem_DAMAS),idx);
fprintf('\n--- DAMAS 0.2 pico: theta=%.1f, phi=%.1f ---\n', theta_vec(tM), phi_vec(pM));
fprintf('--- Caixa A: theta=%.1f, phi=%.1f ---\n', cA_theta, cA_phi);