clear; close all; clc;
%==========================================================================
% Conversor_WPE.m  (roda da PASTA DO PROJETO)
% WPE oficial da NTT (wpe.p, v1.33) nos 6 canais, cortando para DUR_USAR s.
%
% Pre-requisito: em wpe_v1.33/settings/local.m, use blk_len = 5;
%
% Estrutura:
%   IME-RJ_IC_AcousticImaging/
%       Conversor_WPE.m  (rode daqui)
%       gravacoes/mic1.wav ... mic6.wav   (61 s; serao cortados)
%       wpe_v1.33/  (wpe.p, settings/local.m)
%
% Saida: gravacoes/mic1_wpe.wav ... mic6_wpe.wav  (dereverberados, ~15 s)
%==========================================================================

DUR_USAR = 15;   % segundos usados (Rx satura ~10-15 s; suficiente p/ DAMAS)

%--------------------------------------------------------------
% (0) Caminhos
%--------------------------------------------------------------
projroot = fileparts(mfilename('fullpath'));
grav_dir = fullfile(projroot, 'gravacoes');
wpe_dir  = fullfile(projroot, 'wpe_v1.33');
addpath(wpe_dir);

%--------------------------------------------------------------
% (1) Le os 6 canais e corta para DUR_USAR segundos
%--------------------------------------------------------------
arq = {'mic1.wav','mic2.wav','mic3.wav','mic4.wav','mic5.wav','mic6.wav'};
M = numel(arq);

x = [];
for m = 1:M
    [tmp, FS] = audioread(fullfile(grav_dir, arq{m}));
    if m == 1, Ncut = min(length(tmp), round(DUR_USAR*FS)); end
    L = min(Ncut, length(tmp));
    seg = zeros(Ncut,1);
    seg(1:L) = tmp(1:L);
    x = [x, seg];          %#ok<AGROW>  (mics nas colunas)
end
fprintf('Usando %d canais, %d amostras (%.1f s), FS=%d.\n', ...
    M, size(x,1), size(x,1)/FS, FS);

%--------------------------------------------------------------
% (2) Roda o WPE oficial (entra na wpe_v1.33 pois cfgs eh relativo)
%--------------------------------------------------------------
old  = cd(wpe_dir);
cfgs = 'settings/local.m';
disp('Rodando WPE oficial (blocos de 5 s)... aguarde.'); tic
try
    y = wpe(x, cfgs);
catch ME
    cd(old); rethrow(ME);
end
cd(old);
fprintf('WPE concluido em %d s.\n', round(toc));

%--------------------------------------------------------------
% (3) Normaliza e salva
%--------------------------------------------------------------
for m = 1:size(y,2)
    ym = y(:,m); pico = max(abs(ym));
    if pico > 0, ym = 0.98*ym/pico; end
    nome = fullfile(grav_dir, sprintf('mic%d_wpe.wav', m));
    audiowrite(nome, ym, FS);
    fprintf('Salvo: %s\n', nome);
end
disp('Pronto. Ouca um *_wpe.wav e depois aponte o DAMAS para eles.');