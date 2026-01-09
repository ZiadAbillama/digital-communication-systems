#EECE 442: Digital Communication Systems Project

A comprehensive MATLAB implementation of a complete digital communication system, covering signal sampling, quantization, source coding, channel coding, modulation, and equalization techniques.

## 📋 Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Part 1: Source Encoding and Decoding](#part-1-source-encoding-and-decoding)
  - [Signal Representation through Sampling](#signal-representation-through-sampling)
  - [Quantization](#quantization)
  - [Source Coding](#source-coding)
  - [Bonus: Dictionary Coding (LZ77)](#bonus-dictionary-coding-lz77)
- [Part 2: Channel Coding and Communication](#part-2-channel-coding-and-communication)
  - [Channel Coding](#channel-coding)
  - [Modulation and Demodulation](#modulation-and-demodulation)
  - [Channel Simulation](#channel-simulation)
  - [Equalization Techniques](#equalization-techniques)
- [Key Results](#key-results)
- [Requirements](#requirements)
- [Usage](#usage)
- [Mathematical Background](#mathematical-background)
- [References](#references)

## 🎯 Overview

This project implements a full digital communication pipeline from end to end, demonstrating fundamental concepts in signal processing and information theory. The system processes analog signals through sampling, quantization, compression, channel coding, modulation, transmission over noisy channels, equalization, and reconstruction.

### Key Features

- **Nyquist-Shannon sampling theorem verification** with reconstruction
- **Fourier series approximation** with convergence analysis
- **Optimal quantization** using Lloyd-Max algorithm
- **Lossless compression** via Huffman and block Huffman coding
- **Dictionary-based compression** using LZ77
- **Multiple modulation schemes**: BPSK, QPSK, 8-PSK, 16-QAM
- **Channel coding** with repetition codes
- **Advanced equalization**: Viterbi MLSE and Zero-Forcing
- **Comprehensive BER analysis** over AWGN and ISI channels

## 📁 Project Structure

```
442 Project/
├── 442 Part1/                          # Source Encoding/Decoding
│   ├── Signal Representation through Sampling/
│   │   ├── Sampling/
│   │   │   └── codes/
│   │   │       ├── sample.m            # Signal sampling function
│   │   │       └── sampletest.m        # Nyquist sampling tests
│   │   ├── Reconstruct/
│   │   │   └── codes/
│   │   │       ├── reconstruct.m       # Sinc interpolation
│   │   │       └── reconstructtest.m   # Reconstruction validation
│   │   └── FourierSeries/
│   │       └── codes/
│   │           ├── ffs.m               # Finite Fourier Series
│   │           ├── testffs.m           # FFS approximation tests
│   │           └── testffserror.m      # Convergence analysis
│   ├── Quantization/
│   │   ├── Uniform and two level quantization/
│   │   │   └── codes/
│   │   │       ├── quan.m              # Quantization function
│   │   │       ├── quantest.m          # Quantization tests
│   │   │       └── quantestMSE.m       # MSE vs levels analysis
│   │   └── Lloydmax/
│   │       └── codes/
│   │           ├── lloyd_max.m         # Lloyd-Max algorithm
│   │           └── Lloyd_maxtest.m     # Optimal quantizer tests
│   ├── Source coding and decoding/
│   │   ├── Huffman/
│   │   │   └── Codes/
│   │   │       ├── huffman_encode.m    # Huffman encoder
│   │   │       ├── huffman_decode.m    # Huffman decoder
│   │   │       ├── huffman_encode_test.m
│   │   │       └── huffman_decode_test.m
│   │   ├── Blocksourcecoding/
│   │   │   └── codes/
│   │   │       ├── block_huffman.m     # Block Huffman coding
│   │   │       └── block_huffman_testing.m
│   │   └── fullchain/
│   │       └── codes/
│   │           └── fullchain.m         # Complete pipeline integration
│   └── Bonus/
│       └── Codes/
│           ├── lz77_encode.m           # LZ77 encoder
│           └── lz77_encode_test.m      # LZ77 tests
│
└── 442 Part 2/                         # Channel Coding & Communication
    ├── Channel coding and decoding/
    │   └── Codes/
    │       ├── repetition_encode.m     # Repetition encoder
    │       ├── repetition_decode.m     # Majority vote decoder
    │       └── Test_enc_dec.m          # Coding gain analysis
    ├── Modulation and Demodulation/
    │   └── Codes/
    │       ├── bpsk_mod.m              # BPSK modulator
    │       ├── bpsk_demod.m            # BPSK demodulator
    │       ├── qpsk_mod.m              # QPSK modulator (Gray)
    │       ├── qpsk_demod.m            # QPSK demodulator
    │       ├── psk8_mod.m              # 8-PSK modulator
    │       ├── psk8_demod.m            # 8-PSK demodulator
    │       ├── qam16_mod.m             # 16-QAM modulator
    │       ├── qam16_demod.m           # 16-QAM demodulator
    │       └── TestingAll.m            # Modulation scheme tests
    ├── Channel, Noise, Detection/
    │   ├── AWGN Channel/
    │   │   └── Codes/
    │   │       └── AWGN_BER.m          # BER curves over AWGN
    │   └── Viterbi Equalizer (MLSE)/
    │       └── Codes/
    │           └── QPSK_ISI_Viterbi.m  # Viterbi equalizer
    └── BONUS/
        └── Codes/
            └── Zero_forcing_equalizer.m # ZF equalizer comparison
```

## 📊 Part 1: Source Encoding and Decoding

### Signal Representation through Sampling

#### 1.1 Generic Sampling

Implements and verifies the Nyquist-Shannon sampling theorem:

**Theory**: A continuous-time signal x(t) with maximum frequency f_max can be perfectly reconstructed from its samples if sampled at f_s ≥ 2f_max (Nyquist rate).

**Implementation**:
- `sample.m`: Samples continuous signals at specified rates
- `reconstruct.m`: Performs sinc interpolation for signal reconstruction

**Key Results**:
- **Below Nyquist (f_s = 100 Hz)**: Severe aliasing, MSE ≈ 0.5
- **Above Nyquist (f_s = 400 Hz)**: Near-perfect reconstruction, MSE ≈ 10^-15

```matlab
% Example usage
f0 = 100; t = 0:1e-5:0.05;
x = cos(2*pi*f0*t);
fs = 400; % Above Nyquist
[t_sample, x_sample] = sample(t, x, fs);
x_reconstructed = reconstruct(t, x_sample, fs);
```

#### 1.2 Sampling through Fourier Series

Analyzes signal representation using Finite Fourier Series (FFS):

**Theory**: A periodic signal can be approximated as:
```
x(t) ≈ Σ(k=-n to n) c_k * exp(j*2π*k*t/T)
```

**Implementation**: `ffs.m` computes Fourier coefficients and reconstructs signals

**Convergence Analysis**:
- **Varying n (fixed T=2s)**: Error decreases exponentially with more harmonics
- **Varying T (fixed n=50)**: Minimum error at T=0.5s (matches signal period)
- **Gibbs phenomenon**: Observable near discontinuities even with high n

### Quantization

#### 2.1 Uniform Quantization

Maps continuous amplitudes to M discrete levels with equal spacing.

**Theory**: For uniform quantization with step size Δ:
```
MSE_theoretical = Δ²/12
```

**Implementation**: `quan.m` with uniform threshold/level generation

**Results**:
| M (levels) | Measured MSE | Theoretical MSE |
|------------|-------------|-----------------|
| 2          | 0.2500      | 0.2500         |
| 4          | 0.0625      | 0.0625         |
| 8          | 0.0156      | 0.0156         |
| 16         | 0.0039      | 0.0039         |
| 32         | 0.0010      | 0.0010         |

**Key Insight**: MSE decreases by ~6 dB per bit (doubling M halves MSE)

#### 2.2 Lloyd-Max Quantization

Optimal quantizer minimizing MSE for a given source distribution.

**Algorithm**:
1. Initialize levels uniformly
2. Update thresholds: t_k = (a_k + a_{k+1})/2
3. Update levels: a_k = E[X | X ∈ region_k]
4. Iterate until convergence (|Δa| < tolerance)

**Implementation**: `lloyd_max.m` with iterative optimization

**Performance Comparison** (M=8 levels):
- Uniform: MSE = 0.0156
- Lloyd-Max: MSE = 0.0143 (~9% improvement)

**Visual Insight**: Lloyd-Max places levels densely where signal amplitude is most probable (e.g., near zero for sinusoids), reducing distortion.

### Source Coding

#### 3.1 Huffman Coding

Variable-length prefix-free code assigning shorter codes to more frequent symbols.

**Theory**: Average length L̄ satisfies:
```
H(A) ≤ L̄ < H(A) + 1
```
where H(A) = -Σ p_i log₂(p_i) is the entropy.

**Implementation**:
- `huffman_encode.m`: Builds Huffman tree and encodes sequences
- `huffman_decode.m`: Lossless reconstruction from bitstream

**Example Results**:

**Simple Discrete Source** [1 1 1 2 2 3 3 3 3 3]:
- Entropy: 1.4855 bits/symbol
- Huffman: 1.5 bits/symbol (99.03% efficiency)
- Fixed-length: 2 bits/symbol
- Compression gain: 25%

**Uniform Quantized Signal** (M=8):
- Entropy: 2.7271 bits/symbol
- Huffman: 2.8092 bits/symbol (97.08% efficiency)
- Fixed-length: 3 bits/symbol
- Compression gain: 6.36%

**Lloyd-Max Quantized Signal** (M=8):
- Entropy: 2.9635 bits/symbol
- Huffman: 3.0000 bits/symbol (98.78% efficiency)
- Fixed-length: 3 bits/symbol
- Compression gain: 1%

**Key Insight**: Lloyd-Max produces more uniform symbol distribution (higher entropy), leaving less redundancy for Huffman to exploit. Most compression gain comes from optimal quantization itself.

#### 3.2 Block Huffman Coding

Extends Huffman coding to k-symbol blocks, approaching entropy more closely.

**Theory**: As block size k increases:
```
L̄_k/k → H(A)  (efficiency → 100%)
```

**Results**:

| Source Type | k | Entropy/k | L̄/k | Efficiency | |A_k| |
|-------------|---|-----------|------|------------|------|
| Discrete    | 1 | 1.4855    | 1.50 | 99.03%     | 3    |
| Discrete    | 2 | 1.4855    | 1.49 | 99.67%     | 9    |
| Discrete    | 4 | 1.4855    | 1.49 | 100.00%    | 27   |
| Uniform     | 2 | 2.7582    | 2.79 | 98.93%     | 64   |
| Lloyd-Max   | 2 | 2.9410    | 2.96 | 99.33%     | 64   |

**Trade-off**: Higher efficiency requires exponentially larger codebooks (|A_k| = M^k).

### Bonus: Dictionary Coding (LZ77)

Compression through pattern repetition rather than symbol frequency.

**Algorithm**:
1. Maintain sliding window of recently seen data
2. Search for longest match in window
3. Emit (offset, length, next_literal) tokens

**Implementation**: `lz77_encode.m` with configurable window/lookahead sizes

**Results**:

**Test Case 1** - Repetitive sequence [1 2 3] × 10:
- Bits/symbol: 0.967
- Compression ratio: 2.069
- Result: ✓ Successful compression

**Test Case 2** - Random sequence [1 2 3 4 5 6 7 8 9]:
- Bits/symbol: 5.0
- Compression ratio: 0.6
- Result: ✗ Expansion (no patterns to exploit)

**Quantized Signal Application**:

| Quantizer | Huffman (bits/sym) | LZ77 (bits/sym) | LZ77 CR |
|-----------|-------------------|-----------------|---------|
| Uniform   | 2.8212            | 3.8152          | 0.7863  |
| Lloyd-Max | 2.9810            | 4.1908          | 0.7159  |

**Conclusion**: For noisy quantized signals with little temporal correlation, symbol-based coding (Huffman) outperforms pattern-based coding (LZ77). LZ77 excels with repetitive textual or structured data.

## 📡 Part 2: Channel Coding and Communication

### Channel Coding

#### Repetition Code (L=3)

Simple error correction: each bit transmitted L times, decoded by majority vote.

**Implementation**:
- `repetition_encode.m`: Repeats each bit L times
- `repetition_decode.m`: Majority voting (threshold L/2)

**Performance** (SNR = 2 dB, BPSK):
- Uncoded BER: 0.08135
- Coded BER: 0.00640
- Improvement: **12.7× better**

**Trade-off**: Reduces throughput by factor L (from 1 bit/symbol to 1/3 bit/symbol).

### Modulation and Demodulation

Implemented modulation schemes with unit average symbol energy:

#### BPSK (Binary Phase Shift Keying)
- **Mapping**: 0 → -A, 1 → +A
- **Bits/symbol**: 1
- **Energy**: E_s = A²

```matlab
% Example
x = bpsk_mod([0 1], 1);  % [-1.0, +1.0]
b = bpsk_demod(x);       % [0, 1]
```

#### QPSK (Quadrature PSK with Gray Mapping)
- **Mapping**: 
  - [0 0] → (-1-j)/√2
  - [0 1] → (+1-j)/√2
  - [1 1] → (+1+j)/√2
  - [1 0] → (-1+j)/√2
- **Bits/symbol**: 2
- **Energy**: E_s = A²

```matlab
x = qpsk_mod([1 0], 1);  % (-0.707 + 0.707j)
```

#### 8-PSK (8-Phase Shift Keying)
- **Mapping**: 8 equally spaced phases (0°, 45°, ..., 315°)
- **Bits/symbol**: 3
- **Energy**: E_s = A²

#### 16-QAM (16-Quadrature Amplitude Modulation)
- **Mapping**: 4×4 Gray-coded square constellation
- **Normalization**: 1/√10 for unit energy
- **Bits/symbol**: 4
- **Levels**: {-3, -1, +1, +3} × 1/√10

```matlab
x = qam16_mod([0 0 1 1], 1);  % (-0.95 + 0.32j)
```

### Channel Simulation

#### AWGN Channel BER Performance

Bit Error Rate simulation over Additive White Gaussian Noise channel.

**Channel Model**:
```
y = x + n, where n ~ CN(0, N₀)
```

**SNR Definition**: 10log₁₀(E_s/N₀) in dB

**Results Summary**:

| SNR (dB) | BPSK Uncoded | BPSK Coded | QPSK Uncoded | QPSK Coded |
|----------|-------------|------------|--------------|------------|
| 0        | 7.86×10⁻²   | 3.37×10⁻²  | 1.30×10⁻¹    | 7.76×10⁻²  |
| 4        | 1.26×10⁻²   | 3.99×10⁻⁴  | 3.61×10⁻²    | 5.69×10⁻³  |
| 8        | 3.82×10⁻⁴   | 8.00×10⁻⁷  | 2.59×10⁻³    | 2.26×10⁻⁵  |
| 10       | 3.84×10⁻⁵   | 1.00×10⁻⁸  | 4.67×10⁻⁴    | 9.00×10⁻⁷  |

**Key Observations**:
- Repetition coding provides ~1 order of magnitude BER improvement
- BPSK and QPSK have identical BER at same E_b/N₀ (both have same minimum distance)
- 8-PSK and 16-QAM show higher BER (packed constellations → reduced distance)

#### Higher-Order Modulation Comparison

**At 8 dB SNR**:
- BPSK: 3.82×10⁻⁴
- QPSK: 2.59×10⁻³
- 8-PSK: 1.45×10⁻²
- 16-QAM: 1.52×10⁻²

**Trade-off**: Higher spectral efficiency (more bits/symbol) vs. robustness (higher BER).

### Equalization Techniques

#### ISI Channel Model

**Channel**: y_i = 2x_i + x_{i-1} + n_i

This 2-tap channel introduces intersymbol interference (ISI) where the previous symbol affects the current received sample.

**Normalized Channel**: h = [2, 1] / √5 to maintain ||h||² = 1

#### Viterbi MLSE (Maximum Likelihood Sequence Estimation)

Optimal sequence detector using dynamic programming.

**Algorithm**:
1. Define trellis with 4 states (one per QPSK symbol)
2. For each time step, compute branch metrics:
   ```
   metric = |y_k - (h₀*x_k + h₁*x_{k-1})|²
   ```
3. Update path metrics (keep minimum)
4. Traceback to find optimal sequence

**Implementation**: `QPSK_ISI_Viterbi.m`

**Performance**:
- At low SNR (0-4 dB): ISI degrades performance vs. AWGN
- At high SNR (8+ dB): MLSE nearly matches AWGN baseline
- Successfully mitigates ISI through joint detection

#### Zero-Forcing (ZF) Equalizer

Linear filter inverting the channel impulse response.

**Design**:
```
g_ZF = (H^T H)^{-1} H^T e_d
```
where H is the channel convolution matrix and e_d is the desired delayed impulse.

**Implementation**: `Zero_forcing_equalizer.m`

**Performance Comparison** (at 10 dB):
- AWGN baseline: ~10⁻⁵ BER
- Viterbi MLSE: ~10⁻⁴ BER
- Zero-Forcing: ~10⁻³ BER

**Key Insight**: ZF removes ISI perfectly in noiseless conditions but amplifies noise significantly. MLSE outperforms ZF by jointly optimizing over sequences rather than inverting the channel directly.

**Noise Enhancement in ZF**:
When H is ill-conditioned or has spectral nulls, (H^T H)^{-1} amplifies noise components, degrading performance compared to MLSE.

## 🎯 Key Results

### Part 1 Highlights

1. **Nyquist Theorem Verification**: Sampling above 2f_max enables perfect reconstruction (MSE < 10⁻¹⁴)

2. **Quantization Efficiency**: Lloyd-Max achieves 5-10% MSE reduction vs. uniform quantization for non-uniform sources

3. **Compression Performance**:
   - Huffman coding: 97-99% efficiency (L̄ ≈ H)
   - Block Huffman: Approaches 100% efficiency as k increases
   - LZ77: Effective for repetitive data; fails on noisy signals

4. **Full Chain Integration**: 
   - Uniform + Huffman: 2.85 bits/symbol, SNR = 15.68 dB
   - Lloyd-Max + Huffman: 2.98 bits/symbol, SNR = 16.63 dB
   - Trade-off: Lloyd-Max uses 4.6% more bits for 0.95 dB SNR gain

### Part 2 Highlights

1. **Channel Coding Gain**: Repetition code (L=3) provides 12.7× BER improvement at cost of 3× bandwidth

2. **Modulation Trade-offs**:
   - BPSK/QPSK: Most robust (E_b/N₀ advantage)
   - 8-PSK/16-QAM: 3-4× spectral efficiency but 1-2 orders of magnitude higher BER

3. **Equalization Performance**:
   - Viterbi MLSE: Near-optimal for ISI channels
   - Zero-Forcing: 10× worse BER than MLSE due to noise enhancement

4. **System Integration**: Complete pipeline demonstrates fundamental trade-offs:
   - Rate vs. Distortion (quantization)
   - Rate vs. Efficiency (source coding)
   - Rate vs. Reliability (channel coding)
   - Bandwidth vs. Power (modulation)

## 💻 Requirements

- **MATLAB** R2019a or later
- **Required Toolboxes**: None (pure MATLAB implementation)
- **Optional**: Communications System Toolbox (for validation only)

## 🚀 Usage

### Part 1: Source Coding Chain

```matlab
% 1. Sample and reconstruct a signal
f0 = 100; t = 0:1e-5:0.05;
x = cos(2*pi*f0*t);
fs = 400; % Above Nyquist
[t_s, x_s] = sample(t, x, fs);
x_reconstructed = reconstruct(t, x_s, fs);

% 2. Quantize using Lloyd-Max
signal = sin(2*pi*5*t) + 0.3*randn(size(t));
M = 8;
[xq, thr, lvl] = lloyd_max(signal, M);

% 3. Compress with Huffman coding
[encoded_bits, dict, avglen] = huffman_encode(xq);
decoded = huffman_decode(encoded_bits, dict);

% 4. Block Huffman for higher efficiency
k = 2;
[bits_block, dict_block, avglen_block, eff] = block_huffman(xq, k);
```

### Part 2: Communication Chain

```matlab
% 1. Generate random bits
N = 10000;
b = randi([0 1], 1, N);

% 2. Channel coding (repetition)
L = 3;
c = repetition_encode(b, L);

% 3. Modulation (QPSK)
N_sym = length(c)/2;
b_mat = reshape(c, 2, N_sym).';
x = zeros(1, N_sym);
for i = 1:N_sym
    x(i) = qpsk_mod(b_mat(i,:), 1);
end

% 4. AWGN channel
SNR_dB = 8;
SNR_lin = 10^(SNR_dB/10);
sigma = sqrt(1/(2*SNR_lin));
y = x + sigma*(randn(size(x)) + 1j*randn(size(x)));

% 5. Demodulation
c_hat = zeros(size(c));
for i = 1:N_sym
    pair = qpsk_demod(y(i));
    c_hat(2*i-1:2*i) = pair;
end

% 6. Channel decoding
b_hat = repetition_decode(c_hat, L);

% 7. Compute BER
BER = mean(b ~= b_hat);
```

### Running Complete Simulations

```matlab
% Part 1: Full source coding chain
run('442 Part1/Source coding and decoding/fullchain/codes/fullchain.m')

% Part 2: AWGN BER curves
run('442 Part 2/Channel, Noise, Detection/AWGN Channel/Codes/AWGN_BER.m')

% Part 2: ISI channel with Viterbi equalizer
run('442 Part 2/Channel, Noise, Detection/Viterbi Equalizer (MLSE)/Codes/QPSK_ISI_Viterbi.m')
```

## 📐 Mathematical Background

### Sampling Theorem

**Nyquist-Shannon Sampling Theorem**: A bandlimited signal x(t) with maximum frequency f_max can be perfectly reconstructed from samples taken at rate f_s ≥ 2f_max using:

```
x(t) = Σ x[n] · sinc((t - nT_s)/T_s)
```

where T_s = 1/f_s and sinc(x) = sin(πx)/(πx).

### Quantization

**Uniform Quantizer MSE**: For step size Δ and uniform source over [x_min, x_max]:

```
MSE = Δ²/12
Δ = (x_max - x_min)/M
```

**Lloyd-Max Conditions**: Optimal quantizer satisfies:
1. **Nearest-neighbor**: t_k = (a_k + a_{k+1})/2
2. **Centroid**: a_k = E[X | X ∈ region_k]

### Information Theory

**Entropy** (average information content):
```
H(A) = -Σ p_i log₂(p_i)  [bits/symbol]
```

**Huffman Bound**:
```
H(A) ≤ L̄ < H(A) + 1
```

**Block Coding Efficiency**:
```
η_k = (k·H(A))/(L̄_k) → 1 as k → ∞
```

### Channel Coding

**Repetition Code**: Rate R = 1/L

**Majority Vote Decision**: Bit decoded as 1 if more than L/2 symbols are 1.

**Coding Gain** (approximate for high SNR):
```
Gain ≈ 10log₁₀(L) dB
```

### Modulation

**Symbol Energy**: E_s = E[|x|²]

**Bit Energy**: E_b = E_s / log₂(M)

**AWGN BER** (BPSK):
```
BER = Q(√(2E_b/N₀))
```
where Q(x) = (1/√(2π)) ∫_x^∞ exp(-t²/2) dt

### Equalization

**MLSE Metric** (Viterbi):
```
Minimize: Σ |y_k - ĥ*x̂_k|²
```

**Zero-Forcing Solution**:
```
g_ZF = (H^T H)^{-1} H^T e_d
```

Forces ideal impulse response but may amplify noise.

## 📚 References

1. **Sampling & Reconstruction**: 
   - Nyquist, H. (1928). "Certain topics in telegraph transmission theory"
   - Shannon, C. E. (1949). "Communication in the presence of noise"

2. **Quantization**:
   - Lloyd, S. P. (1982). "Least squares quantization in PCM"
   - Max, J. (1960). "Quantizing for minimum distortion"

3. **Source Coding**:
   - Huffman, D. A. (1952). "A method for the construction of minimum-redundancy codes"
   - Ziv, J., & Lempel, A. (1977). "A universal algorithm for sequential data compression"

4. **Digital Communications**:
   - Proakis, J. G., & Salehi, M. (2008). "Digital Communications" (5th ed.)
   - Cover, T. M., & Thomas, J. A. (2006). "Elements of Information Theory" (2nd ed.)

5. **Equalization**:
   - Forney, G. D. (1973). "The Viterbi algorithm"
   - Haykin, S. (2002). "Adaptive Filter Theory" (4th ed.)

## 📄 License

This project was developed for educational purposes as part of EECE 442: Communication Systems.

## 👥 Author

Project developed for EECE 442 - Digital Communication Systems

---

**Note**: All implementations are optimized for clarity and educational value rather than computational efficiency. For production systems, consider using built-in MATLAB Communication Toolbox functions.
