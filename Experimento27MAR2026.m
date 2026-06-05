clear; close all; clc;
%----------------
% Inicializações:
%----------------
c = 343; % velocidade do som em m/s (ajustar conforme temperatura)
f = 700; % frequência do tom num dos alto-falantes
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
P  = [p1; p2; p3; p4; p5; p6]; % Matriz de posições
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

% if size(x1,2) > 1, disp ('Sineal estéreo'); end
% Usar x1 = x1(:,1);
% Usar x2 = x1(:,2);

tam = min([length(x1), length(x2), length(x3), length(x4), length(x5), length(x6)]);
x1 = x1(1:tam); x2 = x2(1:tam); x3 = x3(1:tam);
x4 = x4(1:tam); x5 = x5(1:tam); x6 = x6(1:tam);
t = (0:(tam-1))/fs;
%-----------------------------------
% Sinal x1 no domínio da frequência:
%-----------------------------------
%x1=x1-mean(x1);
% [X1,w]=freqz(x1,1,1000000,'whole');
% plot((fs/2)*(w-pi)/pi/1000,20*log10(abs(fftshift(X1))))
% xlabel('f (kHz)')
% ylabel('|X_1(f)|_{dB}')
% grid
% axis([0 fs/2000 20 100])
% title('Sinal x_1(t) no domínio da frequência')
% return
%-----------------------------------------------------------------------
% Filtrando os sinais por um BPF centrado em 700 Hz
%-----------------------------------------------------------------------
fo=f; 
ORDEM=1500;
deltaf=100;
freqs=[0 fo-deltaf fo-deltaf/2 fo+deltaf/2 fo+deltaf fs/2]/(fs/2);
ganhos=[0 0 1 1 0 0]; 
h = firpm(ORDEM,freqs,ganhos,[1 1 1]);
%[H,w]=freqz(h,1,10000);
%hold on 
%plot((fs/2000)*w/pi,20*log10(abs(H)),'k','linewidth',2)
%plot(freqs*fs/2000,20*log10(abs(ganhos+1e-6)),'r','linewidth',2)
%return
% tic
% x1=filter(h,1,x1);
% toc
% x2=filter(h,1,x2);
% x3=filter(h,1,x3);
% x4=filter(h,1,x4);
% x5=filter(h,1,x5);
% x6=filter(h,1,x6);
x1=filtfilt(h,1,x1); % zero delay
x2=filtfilt(h,1,x2); % zero delay
x3=filtfilt(h,1,x3); % zero delay
x4=filtfilt(h,1,x4); % zero delay
x5=filtfilt(h,1,x5); % zero delay
x6=filtfilt(h,1,x6); % zero delay
%-----------------------------------------------------------------------
% Preparando sinais (tornando-os analíticos) e montando matrix snapshots
%-----------------------------------------------------------------------
% Tornando os sinais analíticos
x1 = hilbert(x1); x2 = hilbert(x2);
[X1,w]=freqz(x1,1,100000,'whole');
plot((fs/2)*(w-pi)/pi/1000,20*log10(abs(fftshift(X1))),'b')
figure
plot((fs/2)*(w-pi)/pi/1000,(abs(fftshift(X1))),'b')
%return
x3 = hilbert(x3); x4 = hilbert(x4);
x5 = hilbert(x5); x6 = hilbert(x6);
snapshots = [x1.'; x2.';x3.'; x4.'; x5.'; x6.']; % M x tam
%return
%-------------------------------------
% Estimando matriz de covariância Rx:
%-------------------------------------
Rx=snapshots*snapshots'/length(snapshots);

% Forward-Backward Averaging
J = fliplr(eye(M));
Rx = (Rx + J * conj(Rx) * J) / 2;

% Diagonal Loading
epsilon = 0.01 * trace(Rx) / M;
Rx = Rx + epsilon * eye(M);

% Inverter Rx
Rxi = inv(Rx);
%return
disp('Be patient and wait!')
tic
% Grid azimute (phi de 0 a 360) e zênite (theta de 0 a 180)
deltagraus=1;
dimensaoy=length(0:deltagraus:180);
dimensaox=length(0:deltagraus:360);
Imagem = zeros(dimensaoy,dimensaox);
contaphi=0;
for phi=0:deltagraus:360
    contaphi=contaphi+1;
    contatheta=0;
    for theta=0:deltagraus:180
        contatheta=contatheta+1;
        % Para cada direção, ache um processador para aquela direção
        %--------
        %  D & S:
        %--------
        % Vetor de direção da fonte:
        vectork = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];
        %  w = zeros(M,1); % Inicializando o vetor do processador
        %  for i = 1:M
        %      % Calculo do tempo de atraso
        %      delta_t = dot(P(i, :), vectork') / c; % Produto escalar
        %      % Cálculo do ganho do processador
        %      w(i) = exp(1i * (2 * pi * f) * delta_t);
        %  end
        % w=w/M
        %w=exp(-1i*2*pi*f*P*vectork/c)/M;
        %pause
        %------
        % MPDR:
        %------
        vaux=exp(-1i*2*pi*f*P*vectork/c);
        w=Rxi*vaux/(vaux'*Rxi*vaux);
        % Obter a saída y(n)=w^H snapshot: y(k) = snapshots.'*conj(w);
        % Medir a variância de y(n) para a dada direção e colocar
        % valor encontrado na matriz Imagem(contaphi,contathete)
        %Imagem(contatheta,contaphi) = var(snapshots.'*conj(w));
        Imagem(contatheta,contaphi) = real(w'*Rx*w); % igual ao anterior
        %---------------------------------------------------------------
     end
end
clc
fprintf('You were patient enough to wait for %d seconds.\n', round(toc));

% Parâmetros
phi = 0:deltagraus:360; % Valores de azimute (0 a 360 graus, de 10 em 10)
theta = 0:deltagraus:180; % Valores de zênite (0 a 180 graus, de 10 em 10)
[PHI, THETA] = meshgrid(phi, theta); % Cria a grade 2D
% Plotando a Imagem
% 2 D:
%figure;
%imagesc(phi, theta, Imagem); % Plota a matriz como uma imagem
% 3 D:
% Heatmap 2D
figure
imagesc(0:deltagraus:360, 0:deltagraus:180, 10*log10(Imagem/max(Imagem(:))))
colorbar; colormap jet
xlabel('\phi (azimute em graus)'); ylabel('\theta (zênite em graus)')
title('MPDR + FB + DL — Pseudoespectro (dB)')
hold on
plot(116.56, 77.60, 'w+', 'MarkerSize', 15, 'LineWidth', 2) % Caixa A
plot(76.02,  83.22, 'wx', 'MarkerSize', 15, 'LineWidth', 2) % Caixa B
legend('Caixa A esperada', 'Caixa B esperada')

% Heatmap 3D
figure
Z = Imagem;
h = surf(PHI,THETA,Z,Imagem,'EdgeColor','none');
axis vis3d; view(45,30); rotate3d on;
xlabel('\phi'); ylabel('\theta')
title('MPDR + FB + DL — Heatmap 3D')

% Pico principal
[maxVal, idx] = max(Imagem(:));
[thetaIdx, phiIdx] = ind2sub(size(Imagem), idx);
fprintf('Pico principal MPDR+FB+DL: phi=%.1f, theta=%.1f\n', (phiIdx-1)*deltagraus, (thetaIdx-1)*deltagraus);