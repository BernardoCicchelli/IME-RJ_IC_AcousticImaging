clear; close all; clc;
%==========================================================================
% DAMAS on-the-fly FINAL — Gauss-Seidel LINHA-A-LINHA (validado correto)
% sobre sinais WPE. Resolve o tempo encolhendo o grid para a regiao de
% interesse (em vez de vetorizar, que quebra a convergencia do DAMAS).
%
% LICAO (p/ relatorio): a vetorizacao por blocos (Jacobi ou GS-por-blocos)
% NAO preserva a convergencia do DAMAS — zera ou distorce o mapa. O loop
% sequencial linha-a-linha (Algoritmo 1 do paper) eh o unico fiel.
%==========================================================================

c=343; f=700; M=6;
p1=(1/100)*[237.0;0.0;105.89].'; p2=(1/100)*[-118.60;205.40;105.89].';
p3=(1/100)*[-118.60;-205.40;105.89].'; p4=(1/100)*[-175.00;0.0;265.89].';
p5=(1/100)*[87.50;-151.60;265.89].'; p6=(1/100)*[87.50;151.60;265.89].';
P=[p1;p2;p3;p4;p5;p6];

%--------------------------------------------------------------
% Sinais WPE + BPF 700 + Hilbert + Rx
%--------------------------------------------------------------
pasta='gravacoes/';
arqs={'mic1_wpe.wav','mic2_wpe.wav','mic3_wpe.wav','mic4_wpe.wav','mic5_wpe.wav','mic6_wpe.wav'};
[x1,fs]=audioread([pasta arqs{1}]); tam=length(x1);
X=zeros(M,tam); X(1,:)=x1(:).';
for m=2:M
    [xm,~]=audioread([pasta arqs{m}]); L=min(tam,length(xm)); X(m,1:L)=xm(1:L).';
end
fo=f; ORDEM=1500; deltaf=100;
freqs=[0 fo-deltaf fo-deltaf/2 fo+deltaf/2 fo+deltaf fs/2]/(fs/2);
h=firpm(ORDEM,freqs,[0 0 1 1 0 0],[1 1 1]);
for m=1:M, xm=filtfilt(h,1,X(m,:).'); X(m,:)=hilbert(xm).'; end
Rx=(X*X')/tam;

%--------------------------------------------------------------
% GRID 0.2 grau — REGIAO DE INTERESSE (em torno das caixas)
% Encolher o grid mantem o loop correto E o tempo baixo.
% Caixa A ~ (78,117). Cobrimos uma janela com folga.
%--------------------------------------------------------------
dg=0.2;
theta_vec=70:dg:90;     % zenite (em torno das caixas)
phi_vec=100:dg:130;     % azimute (em torno da Caixa A)
Ntheta=numel(theta_vec); Nphi=numel(phi_vec); N=Ntheta*Nphi;
fprintf('Grid 0.2: N=%d direcoes (regiao de interesse).\n', N);

%--------------------------------------------------------------
% Asteer + b (DAS)
%--------------------------------------------------------------
disp('Pre-computando...'); tic
Asteer=zeros(M,N); col=0;
for it=1:Ntheta
    th=theta_vec(it); st=sind(th); ct=cosd(th);
    for ip=1:Nphi
        ph=phi_vec(ip); col=col+1;
        Asteer(:,col)=exp(-1i*2*pi*f*P*[st*cosd(ph);st*sind(ph);ct]/c);
    end
end
RA=Rx*Asteer;
b_img=real(sum(conj(Asteer).*RA,1)).'/(M^2);
fprintf('Pronto em %d s.\n', round(toc));

%--------------------------------------------------------------
% DAMAS on-the-fly — LINHA-A-LINHA (Gauss-Seidel puro, fiel ao Algoritmo 1)
% A i-esima linha de A eh gerada on-the-fly; A nunca eh montada inteira.
%--------------------------------------------------------------
NLOOPS=100;
p=zeros(N,1);
% diagonal: |a_i^H a_i|^2/M^2 = M^2/M^2 = 1 para steering nao-normalizado,
% mas calculo explicito p/ robustez:
normcol = real(sum(conj(Asteer).*Asteer,1)).';   % = M para cada coluna
Adiag = (normcol.^2)/(M^2);                       % = 1

fprintf('Rodando DAMAS linha-a-linha (%d iter, N=%d)...\n', NLOOPS, N);
tic
for loop=1:NLOOPS
    for i=1:N
        g = Asteer(:,i)' * Asteer;        % 1 x N (linha i, on-the-fly)
        Ai = (abs(g).^2)/(M^2);
        interf = Ai*p - Ai(i)*p(i);
        p(i) = max(0, (b_img(i) - interf)/Adiag(i));
    end
    if mod(loop,10)==0, fprintf('  iter %d/%d (%.0f s)\n',loop,NLOOPS,toc); end
end
fprintf('DAMAS concluido em %d s (%.1f min).\n', round(toc), toc/60);

Imagem_DAS=reshape(b_img,Nphi,Ntheta).';
Imagem_DAMAS=reshape(p,Nphi,Ntheta).';

%--------------------------------------------------------------
% Plots
%--------------------------------------------------------------
cA_phi=116.56; cA_theta=77.60; cB_phi=76.50; cB_theta=83.22;
figure
imagesc(phi_vec,theta_vec,Imagem_DAS/max(Imagem_DAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi'); ylabel('\theta'); title('DAS 0.2° — WPE (zoom)')
hold on; plot(cA_phi,cA_theta,'w+','MarkerSize',15,'LineWidth',2)

figure
imagesc(phi_vec,theta_vec,Imagem_DAMAS/max(Imagem_DAMAS(:)))
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi'); ylabel('\theta'); title('DAMAS on-the-fly 0.2° — WPE')
hold on; plot(cA_phi,cA_theta,'w+','MarkerSize',15,'LineWidth',2)

[~,idx]=max(Imagem_DAMAS(:)); [tM,pM]=ind2sub(size(Imagem_DAMAS),idx);
fprintf('\n--- DAMAS 0.2 pico: theta=%.1f, phi=%.1f ---\n', theta_vec(tM), phi_vec(pM));
fprintf('--- Caixa A: theta=%.1f, phi=%.1f ---\n', cA_theta, cA_phi);
save('damas_02_resultado.mat','Imagem_DAMAS','Imagem_DAS','theta_vec','phi_vec');