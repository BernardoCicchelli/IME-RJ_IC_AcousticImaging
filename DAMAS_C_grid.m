clear; close all; clc;
%==========================================================================
% DAMAS-C em DOIS ESTAGIOS — localizacao refinada das fontes coerentes
% (Brooks & Humphreys, AIAA 2006-2654), array 3D de 6 mics, projeto IME.
%
% Estagio 1 (GROSSO): zona A larga em torno das DUAS caixas, passo dg1.
%                     Acha os 2 picos (sem usar as coords esperadas como
%                     centro fino -> resultado honesto).
% Estagio 2 (FINO):   roda DAMAS-C de novo numa janela ESTREITA em torno
%                     de CADA pico do estagio 1, com passo dg2 pequeno.
% Pos-processo:       interpolacao parabolica 2D (sub-grid) p/ coordenada
%                     continua de cada fonte.
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

%---------------- le os 6 canais WPE ----------------
pasta = 'gravacoes/';
arqs = {'mic1_wpe.wav','mic2_wpe.wav','mic3_wpe.wav', ...
        'mic4_wpe.wav','mic5_wpe.wav','mic6_wpe.wav'};
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
G = (X*X')/tam;                       % CSM 6x6

cA_theta=77.60; cA_phi=116.56;  cB_theta=83.22; cB_phi=76.50;

%==========================================================================
% FUNCAO LOCAL (definida no fim) damasc_zone() roda DAMAS-C numa lista de
% pontos (theta,phi) e devolve auto-potencias + matriz Xc completa.
% Aqui chamamos ela duas vezes.
%==========================================================================

%--------------------------------------------------------------------------
% ESTAGIO 1 — GROSSO
%--------------------------------------------------------------------------
dg1 = 2; jan1 = 8; NLOOP1 = 300;
[th1, ph1] = make_zone([cA_theta cB_theta], [cA_phi cB_phi], jan1, dg1, ...
                       [60 120], [40 140]);
fprintf('Estagio 1 (grosso): %d pontos.\n', numel(th1));
[auto1, ~] = damasc_zone(th1, ph1, P, G, M, f, c, NLOOP1);

[~, ord] = sort(auto1,'descend');
pk_th = [th1(ord(1)) th1(ord(2))];
pk_ph = [ph1(ord(1)) ph1(ord(2))];
fprintf('  pico 1: theta=%.0f phi=%.0f | pico 2: theta=%.0f phi=%.0f\n', ...
        pk_th(1),pk_ph(1),pk_th(2),pk_ph(2));

%--------------------------------------------------------------------------
% ESTAGIO 2 — FINO em torno de CADA pico (janela estreita, passo pequeno)
%--------------------------------------------------------------------------
dg2 = 1; jan2 = 2; NLOOP2 = 300;   % leve: ~5x5 pts/fonte, AC ~625x625
res = zeros(2,2);    % [theta phi] refinado de cada fonte
for s = 1:2
    [thf, phf] = make_zone(pk_th(s), pk_ph(s), jan2, dg2, [60 120],[40 140]);
    [autof, ~] = damasc_zone(thf, phf, P, G, M, f, c, NLOOP2);
    % interpolacao parabolica 2D sub-grid em torno do max
    [tc, pc] = subgrid_peak(thf, phf, autof, dg2);
    res(s,:) = [tc pc];
    fprintf('Estagio 2 fonte %d: theta=%.2f phi=%.2f (%d pts)\n', ...
            s, tc, pc, numel(thf));
end

%==========================================================================
% RESULTADO FINAL — associa cada fonte refinada a caixa mais proxima
%==========================================================================
exp_th = [cA_theta cB_theta]; exp_ph = [cA_phi cB_phi]; nome={'A','B'};
fprintf('\n================ RESULTADO REFINADO ================\n');
usado = [false false];
for s = 1:2
    d = (res(s,1)-exp_th).^2 + (res(s,2)-exp_ph).^2;
    d(usado) = inf; [~,k]=min(d); usado(k)=true;
    fprintf('Fonte %d -> Caixa %s: theta=%.2f (esp %.1f) | phi=%.2f (esp %.1f)\n',...
        s, nome{k}, res(s,1), exp_th(k), res(s,2), exp_ph(k));
end

%==========================================================================
% FUNCOES LOCAIS
%==========================================================================
function [th, ph] = make_zone(cth, cph, jan, dg, thlim, phlim)
% gera lista de pontos (theta,phi) dentro de +-jan em torno de cada centro
    tv = thlim(1):dg:thlim(2);
    pv = phlim(1):dg:phlim(2);
    [TH,PH] = meshgrid(tv,pv); TH=TH(:); PH=PH(:);
    mask = false(size(TH));
    for k=1:numel(cth)
        mask = mask | (abs(TH-cth(k))<=jan & abs(PH-cph(k))<=jan);
    end
    th = TH(mask); ph = PH(mask);
end

function [auto, Xmat] = damasc_zone(th, ph, P, G, M, f, c, NLOOPS)
% roda DAMAS-C (forma fatorada via kron) numa lista de pontos
    Na = numel(th);
    E = zeros(M, Na);
    for n=1:Na
        vk=[sind(th(n))*cosd(ph(n)); sind(th(n))*sind(ph(n)); cosd(th(n))];
        E(:,n)=exp(-1i*2*pi*f*P*vk/c);
    end
    M2 = M^2;
    Gmat = E'*E;
    Yc = (E'*G*E)/M2;  yC = Yc(:);
    AC = kron(Gmat, Gmat.')/M2;
    Ncc = Na*Na;
    Xc = zeros(Ncc,1); ACd = diag(AC);
    for loop=1:NLOOPS
        for k=1:Ncc
            interf = AC(k,:)*Xc - AC(k,k)*Xc(k);
            r = real((yC(k)-interf)/ACd(k));
            if r<0, r=0; end
            Xc(k)=r;
        end
    end
    Xmat = reshape(Xc,Na,Na);
    auto = real(diag(Xmat));
end

function [tc, pc] = subgrid_peak(th, ph, val, dg)
% interpolacao parabolica 2D em torno do ponto de maximo (sub-grid)
    [~,i] = max(val);
    t0 = th(i); p0 = ph(i);
    % vizinhos ao longo de theta (mesmo phi)
    function v = getv(tt,pp)
        j = find(abs(th-tt)<1e-6 & abs(ph-pp)<1e-6,1);
        if isempty(j), v=NaN; else v=val(j); end
    end
    % theta
    vm=getv(t0-dg,p0); v0=getv(t0,p0); vp=getv(t0+dg,p0);
    dt=0;
    if ~isnan(vm)&&~isnan(vp)
        den=(vm-2*v0+vp); if abs(den)>eps, dt=0.5*(vm-vp)/den; end
    end
    % phi
    wm=getv(t0,p0-dg); w0=v0; wp=getv(t0,p0+dg);
    dp=0;
    if ~isnan(wm)&&~isnan(wp)
        den=(wm-2*w0+wp); if abs(den)>eps, dp=0.5*(wm-wp)/den; end
    end
    tc = t0 + dt*dg;
    pc = p0 + dp*dg;
end