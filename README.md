# 🎙️ Acoustic Imaging — IME-RJ

Projeto de Iniciação Científica desenvolvido no **Instituto Militar de Engenharia (IME-RJ)** na área de **Processamento Digital de Sinais (PDS)**.

---

## 📌 Objetivo

Localizar visualmente fontes sonoras utilizando um array de microfones, gerando um **mapa de calor (heatmap)** que indica a direção de onde o som está vindo — técnica conhecida como **Acoustic Imaging**.

---

## 🔧 Setup Experimental

### Array de Microfones
6 microfones dispostos em duas tríades triangulares em alturas diferentes:

| Microfone | X (m) | Y (m) | Z (m) |
|-----------|-------|-------|-------|
| Mic 1 | 2.370 | 0.000 | 1.059 |
| Mic 2 | -1.186 | 2.054 | 1.059 |
| Mic 3 | -1.186 | -2.054 | 1.059 |
| Mic 4 | -1.750 | 0.000 | 2.659 |
| Mic 5 | 0.875 | -1.516 | 2.659 |
| Mic 6 | 0.875 | 1.516 | 2.659 |

### Fontes Sonoras
Duas caixas de som emitindo tom puro de **700 Hz**:

| Fonte | X (cm) | Y (cm) | Z (cm) | φ (azimute) | θ (zênite) |
|-------|--------|--------|--------|-------------|------------|
| Caixa A | -61.0 | 122.0 | 30.0 | 116.56° | 11.86° |
| Caixa B | 61.0 | 245.0 | 30.0 | 76.50° | 6.27° |

### Câmera
- Lente posicionada na origem (0, 0, 0)
- Modelo: Samsung Galaxy A71
- Os ângulos φ e θ são calculados em relação à lente da câmera

### Ambiente
Laboratório em ambiente **não anecoico** — sem tratamento acústico, com reflexões de paredes, piso e teto.

---

## 🧠 Algoritmos Implementados

### 1. MPDR (Minimum Power Distortionless Response)
Também conhecido como Capon Beamformer. Minimiza a potência total recebida mantendo ganho unitário na direção de interesse.

**Limitação encontrada:** O MPDR foi projetado para uma fonte dominante. Com duas fontes tocando a mesma frequência simultaneamente (fontes correlacionadas) em ambiente não anecoico, o algoritmo gerou múltiplos picos falsos sem identificar corretamente as fontes.

**Referência:** Capon, J. (1969). *High-resolution frequency-wavenumber spectrum analysis*. Proceedings of the IEEE.

---

### 2. MUSIC + Forward-Backward Averaging
O MUSIC (MUltiple SIgnal Classification) é um algoritmo de alta resolução que separa o espaço de sinais do espaço de ruído através da decomposição em autovalores da matriz de covariância.

**Por que Forward-Backward Averaging?**
Como as duas fontes emitem a mesma frequência (700 Hz), os sinais são perfeitamente correlacionados — isso colapsa a matriz de covariância para apenas 1 autovalor significativo, quando deveríamos ter 2. O Forward-Backward Averaging resolve isso:

```
Rx_fb = (Rx + J * conj(Rx) * J) / 2
```

**Resultado obtido:**
- 2 autovalores significativos identificados corretamente ✅
- Pico principal encontrado em φ=106°, θ=16°
- Erro de ~10° em azimute e ~4° em zênite em relação à Caixa A
- Caixa B não localizada com precisão devido à proximidade angular das fontes (40° de separação)

**Referência:** Schmidt, R.O. (1986). *Multiple emitter location and signal parameter estimation*. IEEE Transactions on Antennas and Propagation.

---

## 🔄 Pipeline de Processamento

```
Gravações .wav (6 microfones)
    ↓
Filtro passa-banda FIR (ordem 1500, centrado em 700 Hz, zero delay)
    ↓
Transformada de Hilbert (sinal analítico)
    ↓
Matriz de snapshots (6 × N)
    ↓
Matriz de covariância Rx
    ↓
Forward-Backward Averaging
    ↓
Decomposição em autovalores (MUSIC)
    ↓
Pseudoespectro MUSIC (grid φ × θ)
    ↓
Heatmap 2D e 3D
```

---

## 📁 Estrutura do Repositório

```
IME-RJ_IC_AcousticImaging/
├── Experimento27MAR2026.m      # Código original com MPDR
├── Experimento_MUSIC.m         # Código com MUSIC + Forward-Backward Averaging
└── README.md
```

> ⚠️ Os arquivos de gravação `.wav` não estão incluídos no repositório por serem dados experimentais privados. Para rodar o código, coloque os arquivos `mic1.wav` até `mic6.wav` dentro de uma pasta chamada `gravacoes/`.

---

## 🛠️ Requisitos

- MATLAB R2024a ou superior
- Signal Processing Toolbox
- DSP System Toolbox

---

## 📊 Resultados

### MPDR
Heatmap com múltiplos picos falsos espalhados — algoritmo não conseguiu localizar as fontes no cenário de duas fontes correlacionadas em ambiente não anecoico.

### MUSIC + Forward-Backward Averaging
![MUSIC 2D](resultados/music_2d.png)
![MUSIC 3D](resultados/music_3d.png)

Pico principal encontrado próximo à Caixa A com erro de ~10° em azimute.

---

## 🔮 Próximos Passos

- Regravar as caixas separadamente para eliminar a correlação entre fontes
- Sobrepor o heatmap na imagem fotográfica do ambiente

---

## 👨‍🎓 Autores

Bernardo Cicchelli e Patrick — Alunos de Iniciação Científica  
Orientador: Professor Apolinário  
Instituto Militar de Engenharia (IME-RJ)
