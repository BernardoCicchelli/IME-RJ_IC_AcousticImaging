clear; close all; clc;
%==========================================================================
% DAMAS-C — Deconvolution Approach for the Mapping of Acoustic Sources
%           with COHERENCE (Brooks & Humphreys, AIAA 2006-2654)
%
% Adaptacao do paper (array planar aeroacustico) para o array 3D de 6 mics
% do projeto IME, com DUAS fontes coerentes (Caixas A e B) em 700 Hz.
%
% IDEIA (vs DAMAS comum):
%   DAMAS:    b = A p      -> N incognitas REAIS  (so auto-potencias)
%   DAMAS-C:  YC = AC XC   -> incognitas COMPLEXAS incluindo cross-spectra
%                            X(n0,n) entre pares de pontos => captura coerencia
%
% Ingredientes do paper usados aqui:
%   - Cross-beamforming  Y(n0,n) = e_n0' * G * e_n / M^2        (Eq. 5)
%   - Bloco []_{n0,n} = e_n0 * e_n'  (outer product)            (Eq. 14)
%   - AC(n0n, n0'n') = e_n0' * []_{n0'n'} * e_n / M^2           (Eq. 19)
%   - Iteracao Gauss-Seidel complexa                            (Eq. 36)
%   - Restricao IN-PHASE: Im(X)=0, Re(X)>=0 a cada passo   (sec. Positivity)
%
% ZONING (essencial p/ caber na memoria): em vez do grid inteiro, define-se
% uma ZONA A pequena de candidatos coerentes em torno das duas caixas.
% Custo ~ Na^2 (Na = pontos da zona A). Com Na ~ 60-120 roda em laptop.
%==========================================================================

c = 343; f = 700; M = 6;

%---------------- posicoes dos mics (m) ----------------
p1 = (1/100)*[237.0;   0.0;    105.89].';
p2 = (1/100)*[-118.60; 205.40; 105.89].';
p3 = (1/100)*[-118.60; -205.40;105.89].';
p4 = (1/100)*[-175.00; 0.0;    265.89].';
p5 = (1/100)*[87.50;  -151.60; 265.89].';
p6 = (1/100)*[87.50;   151.60; 265.89].';
P  = [p1; p2; p3; p4; p5; p6];

%---------------- le os 6 canais (use WPE p/ tirar reverb) ----------------
pasta = 'gravacoes/';
arqs = {'mic1_wpe.wav','mic2_wpe.wav','mic3_wpe.wav', ...
        'mic4_wpe.wav','mic5_wpe.wav','mic6_wpe.wav'};
% se nao tiver WPE, troque pelos mic1.wav ... mic6.wav
[x1, fs] = audioread([pasta arqs{1}]);
tam = length(x1);
X = zeros(M, tam); X(1,:) = x1(:).';
for m = 2:M
    [xm, fsm] = audioread([pasta arqs{m}]);
    if fsm ~= fs, error('fs diferente.'); end
    L = min(tam, length(xm)); X(m,1:L) = xm(1:L).';
end

%---------------- BPF 700 Hz + analitico ----------------
fo=f; ORDEM=1500; deltaf=100;
fr=[0 fo-deltaf fo-deltaf/2 fo+deltaf/2 fo+deltaf fs/2]/(fs/2);
h=firpm(ORDEM,fr,[0 0 1 1 0 0],[1 1 1]);
for m=1:M, xm=filtfilt(h,1,X(m,:).'); X(m,:)=hilbert(xm).'; end

% CSM (G no paper) 6x6
G = (X*X')/tam;

%==========================================================================
% GRID + ZONING
%   Grid de varredura amplo so p/ visualizar (DAS). A ZONA A (coerente) e
%   um sub-conjunto pequeno em torno das duas caixas, onde o DAMAS-C roda.
%==========================================================================
dg = 2;                          % passo do grid (2 graus -> mais leve)
theta_vec = 60:dg:120;
phi_vec   = 40:dg:140;
[TH,PH] = meshgrid(theta_vec, phi_vec);
TH = TH(:); PH = PH(:);          % lista de pontos do grid
Ngrid = numel(TH);

% steering de cada ponto do grid (3D): e_n (Mx1)
Eall = zeros(M, Ngrid);
for n = 1:Ngrid
    vk = [sind(TH(n))*cosd(PH(n)); sind(TH(n))*sind(PH(n)); cosd(TH(n))];
    Eall(:,n) = exp(-1i*2*pi*f * P * vk / c);
end

% --- DAS comum (auto-beamform) so p/ referencia visual ---
das = zeros(Ngrid,1);
for n=1:Ngrid, w=Eall(:,n)/M; das(n)=real(w'*G*w); end

%--------------------------------------------------------------------------
% ZONA A: pega os pontos do grid dentro de uma janela em torno de CADA caixa
% (coerencia permitida SO aqui). Mantem Na pequeno!
%--------------------------------------------------------------------------
cA_theta=77.60; cA_phi=116.56;  cB_theta=83.22; cB_phi=76.50;
janela = 8;   % +-8 graus em torno de cada caixa
inA = ( abs(TH-cA_theta)<=janela & abs(PH-cA_phi)<=janela ) | ...
      ( abs(TH-cB_theta)<=janela & abs(PH-cB_phi)<=janela );
idxA = find(inA);            % indices dos pontos da zona A
Na = numel(idxA);
E = Eall(:, idxA);          % steering da zona A (M x Na)
thA = TH(idxA); phA = PH(idxA);
fprintf('Zona A: %d pontos (de %d no grid). AC tera %d x %d.\n', ...
        Na, Ngrid, Na*Na, Na*Na);

%==========================================================================
% MONTA YC e AC  (forma reduzida: indices de pares (a,b), a,b em 1..Na)
%   Yc(a,b) = e_a' G e_b / M^2                         (cross-beamform, Eq.5)
%   Xc(a,b) = cross-spectrum entre fontes a e b (incognita)
%   AC[(a,b),(a',b')] = (e_a' (e_a' ' ... )) -> ver Eq.14/19:
%       []_{a'b'} = e_{a'} e_{b'}'
%       AC = ( e_a' * []_{a'b'} * e_b ) / M^2
%                = (e_a' e_{a'})(e_{b'}' e_b)/M^2
%   => fatoravel! AC[(a,b),(a',b')] = (E'E)(a,a') * (E'E)(b',b) / M^2
%      onde Gmat = E'E (Na x Na). Isso evita o laco quadruplo.
%==========================================================================
Gmat = E' * E;                       % Na x Na  (e_i' e_j)
M2 = M^2;

% Yc (Na x Na) complexo
Yc = (E' * G * E) / M2;              % e_a' G e_b / M^2

% Vetoriza: indice linear k = (a-1)*Na + b  -> par (a,b)
% AC e (Na^2 x Na^2):  AC(k, k') = Gmat(a,a') * conj? cuidado com transpostos.
% Do paper: []_{a'b'} = e_{a'} e_{b'}^{1*}... usamos a forma fatorada:
%   AC[(a,b),(a',b')] = (e_a' e_{a'}) * (e_{b'}' e_b) / M2
%                      = Gmat(a,a') * Gmat(b',b) / M2
% Montagem via kron (Gmat (x) Gmat^T):
AC = kron(Gmat, Gmat.') / M2;        % (Na^2 x Na^2) complexo

yC = Yc(:);                          % (Na^2 x 1)  -- ordem: a varia rapido? 
% MATLAB (:) empilha por COLUNA: k=(b-1)*Na + a  => par (a,b) com a rapido.
% kron(Gmat,Gmat.') segue a mesma convencao, entao yC e AC sao consistentes.

%==========================================================================
% ITERACAO DAMAS-C (Gauss-Seidel complexo) com restricao IN-PHASE
%   X(k) <- ( Yc(k) - sum_{j~=k} AC(k,j) X(j) ) / AC(k,k)
%   AC(k,k) = 1 (Eq.23), mas dividimos pelo valor real p/ robustez.
%   Restricao: para termos diagonais (a==b): Im=0, Re>=0.
%             para cross (a~=b) in-phase: Im=0, Re>=0.
%==========================================================================
Ncc = Na*Na;
Xc = zeros(Ncc,1);
ACdiag = diag(AC);
% mapa k -> (a,b)
[bb, aa] = ind2sub([Na Na], (1:Ncc)');   % por causa do (:) por coluna: a rapido
% (aa = linha = primeiro indice; bb = coluna = segundo indice)

NLOOPS = 300;
fprintf('Rodando DAMAS-C (%d iter, %d incognitas)...\n', NLOOPS, Ncc); tic
for loop = 1:NLOOPS
    for k = 1:Ncc
        interf = AC(k,:)*Xc - AC(k,k)*Xc(k);
        val = (yC(k) - interf) / ACdiag(k);
        % restricao in-phase: parte real >=0, imaginaria zerada
        r = real(val);
        if r < 0, r = 0; end
        Xc(k) = r;          % Im=0 imposto
    end
end
fprintf('DAMAS-C pronto em %d s.\n', round(toc));

% Auto-potencias (a==b) = forca da fonte em cada ponto da zona A
Xmat = reshape(Xc, Na, Na);
auto = real(diag(Xmat));          % Na x 1  -> "mapa" de fontes

%==========================================================================
% RESULTADOS / PLOTS
%==========================================================================
% mapa DAS no grid completo
DAS_img = reshape(das, numel(phi_vec), numel(theta_vec));
figure
imagesc(phi_vec, theta_vec, (DAS_img/max(DAS_img(:))).')
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi'); ylabel('\theta'); title('DAS (referencia)')
hold on
plot(cA_phi,cA_theta,'w+','MarkerSize',15,'LineWidth',2)
plot(cB_phi,cB_theta,'wx','MarkerSize',15,'LineWidth',2)

% mapa DAMAS-C: scatter dos pontos da zona A coloridos pela auto-potencia
figure
scatter(phA, thA, 120, auto/max(auto+eps), 'filled')
set(gca,'YDir','normal'); colorbar; colormap jet
xlabel('\phi'); ylabel('\theta')
title('DAMAS-C — auto-potencias na Zona A (sinais WPE)')
hold on
plot(cA_phi,cA_theta,'k+','MarkerSize',15,'LineWidth',2)
plot(cB_phi,cB_theta,'ko','MarkerSize',15,'LineWidth',2)
xlim([40 140]); ylim([60 120]); grid on
legend('fontes DAMAS-C','Caixa A','Caixa B')

%==========================================================================
% PICOS: dois maiores da auto-potencia
%==========================================================================
[as, ord] = sort(auto, 'descend');
fprintf('\n--- DAMAS-C: 2 maiores fontes na Zona A ---\n');
for k = 1:min(2,Na)
    j = ord(k);
    fprintf('  #%d: theta=%.0f, phi=%.0f   (nivel rel=%.2f)\n', ...
            k, thA(j), phA(j), as(k)/max(as+eps));
end
fprintf('--- Caixa A esperada: theta=%.1f, phi=%.1f ---\n', cA_theta, cA_phi);
fprintf('--- Caixa B esperada: theta=%.1f, phi=%.1f ---\n', cB_theta, cB_phi);

% coerencia recuperada entre os dois picos (diagnostico)
if Na>=2
    a1=ord(1); a2=ord(2);
    Xab = Xmat(a1,a2);
    gamma2 = abs(Xab)^2 / (real(Xmat(a1,a1))*real(Xmat(a2,a2)) + eps);
    fprintf('\nCoerencia^2 estimada entre os 2 picos: %.3f\n', gamma2);
    fprintf('(perto de 1 => DAMAS-C confirma que sao coerentes)\n');
end