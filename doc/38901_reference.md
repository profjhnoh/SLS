# 3GPP TR 38.901 Quick Reference for Tbps_DLSLS Simulator

> **Source:** 3GPP TR 38.901 V17.1.0 (2024-09) — Study on channel model for frequencies from 0.5 to 100 GHz
> **Purpose:** 시뮬레이터 개발 시 빠른 참조용. 전체 원문은 `doc/38901-j10.docx` 참고.

---

## 목차

1. [코드-Spec 매핑](#1-코드-spec-매핑)
2. [채널 생성 절차 (Section 7.5)](#2-채널-생성-절차-section-75)
3. [Path Loss 모델 (Section 7.4.1)](#3-path-loss-모델-section-741)
4. [LOS 확률 모델 (Section 7.4.2)](#4-los-확률-모델-section-742)
5. [O2I Penetration Loss (Section 7.4.3)](#5-o2i-penetration-loss-section-743)
6. [안테나 모델링 (Section 7.3)](#6-안테나-모델링-section-73)
7. [LSP 파라미터 테이블 (Table 7.5-6)](#7-lsp-파라미터-테이블-table-756)
8. [ZSD/ZOD 오프셋 (Table 7.5-7~12)](#8-zsdzod-오프셋-table-757~12)
9. [Scaling Factor 테이블](#9-scaling-factor-테이블)
10. [Calibration 파라미터 (Section 7.8)](#10-calibration-파라미터-section-78)

---

## 1. 코드-Spec 매핑

시뮬레이터 소스 파일과 38.901 섹션의 대응 관계:

| 소스 파일 | 주요 함수/역할 | 38.901 Section |
|---|---|---|
| `src/main.cpp` | 전체 시뮬레이션 루프 | 7.5 전체 흐름 |
| `src/Initiallization.cpp` | BS/MS 배치, 안테나 좌표 계산 | Step 1 (7.5), 7.2, 7.3 |
| `src/channel.cpp` | 채널 계수 생성, LOS/NLOS 판정 | Step 2,11 (7.5), 7.4.2 |
| `src/Link.cpp` | Path loss 계산, O2I loss, 안테나 gain | 7.4.1, 7.4.3 |
| `src/generateLSP.cpp` | DS, ASD, ASA, ZSD, ZSA, K, SF 생성 | Step 4 (7.5), Table 7.5-6 |
| `src/generateSSP.cpp` | 클러스터 delay, power, 각도 생성 | Step 5~9 (7.5) |
| `src/channel_update.cpp` | 채널 계수 시간 업데이트 | Step 11 Doppler |
| `src/setSimulParam.cpp` | cfg 파라미터 파싱 | - |
| `src/setFastfadingParam.cpp` | 상관 행렬 설정 | Step 4 correlation |
| `src/scheduling.cpp` | MU-MIMO 스케줄링 | - |
| `src/receive_downlink.cpp` | 수신 처리, SINR/BLER | - |
| `src/Logging_Point.cpp` | CDF 통계 수집 | 7.8 calibration metrics |

---

## 2. 채널 생성 절차 (Section 7.5)

### Step 1: 환경/레이아웃/안테나 설정
→ `Initiallization.cpp`: `InitializeSystem()`, `Initialdrop()`

- 시나리오 선택 (UMa, UMi, InH, RMa, InF)
- BS/UT 3D 위치 결정, LOS 각도 계산
- 안테나 field pattern `F_rx`, `F_tx` 설정 (GCS 기준)
- BS orientation: Ω_BS,α (bearing), Ω_BS,β (downtilt), Ω_BS,γ (slant)
- UT orientation: Ω_UT,α (bearing), Ω_UT,β (downtilt), Ω_UT,γ (slant)

### Step 2: 전파 조건 할당
→ `channel.cpp`: LOS/NLOS 판정

- Table 7.4.2-1에 따라 LOS/NLOS 결정 (BS-UT 링크 간 uncorrelated)
- Indoor/outdoor 상태 할당 (한 UT의 모든 링크 동일)

### Step 3: Path Loss 계산
→ `Link.cpp`: `Get_CouplingLoss()`

- Table 7.4.1-1 수식으로 각 BS-UT 링크의 path loss 계산

### Step 4: Large Scale Parameter 생성
→ `generateLSP.cpp`: `Generate_LSP()`

- DS, ASD, ASA, ZSD, ZSA, K, SF 생성
- Cross-correlation: Table 7.5-6의 상관 행렬 → Cholesky 분해
- LSP 벡터 순서: s_M = [s_SF, s_K, s_DS, s_ASD, s_ASA, s_ZSD, s_ZSA]^T
- Co-sited sector → 동일 UT에 대해 동일 LSP
- 제한: ASA, ASD ≤ 104°, ZSA, ZSD ≤ 52°

### Step 5: 클러스터 Delay 생성
→ `generateSSP.cpp`

지수 분포에서 delay 추출:
```
τ_n = -r_τ · DS · ln(X_n),  X_n ~ Uniform(0,1)     ... (7.5-1)
τ'_n = sort(τ_n - min(τ))                            ... (7.5-2)
```

LOS인 경우 추가 스케일링:
```
c_τ = sqrt(K_R / (K_R + 1))                          ... (7.5-3)  (하지만 Eq. 7.5-1의 1/λ_c에 주의)
τ''_n = c_τ · τ'_n  (클러스터 power 생성에는 사용 안 함!)   ... (7.5-4)
```

### Step 6: 클러스터 Power 생성

```
P_n = exp(-(τ'_n / (r_τ · c_DS)) - ζ_n)              ... (7.5-5)  (ζ_n: per-cluster shadowing [dB])
P_n = P_n / Σ P_i                                    ... (7.5-6)  (정규화)
```

LOS인 경우:
```
P_LOS = K_R / (K_R + 1)                              ... (7.5-7)
P_n = (1 - P_LOS) · P_n + P_LOS · δ(n-1)            ... (7.5-8)
```

- 각 ray의 power = P_n / M (M = rays per cluster)
- 최대 대비 -25 dB 미만 클러스터 제거

### Step 7: 도래각/출발각 생성

**AOA** (wrapped Gaussian):
```
φ_n,AOA = (2/C_φ) · sqrt(-ln(P_n)/π) · ASA           ... (7.5-9)
```

C_φ^NLOS 값: Table 7.5-2 참조. LOS인 경우 K-factor 보정:
```
C_φ = C_φ^NLOS · (1.1035 - 0.028K - 0.002K²)         ... (7.5-10)
```

부호 및 랜덤 오프셋 추가:
```
φ_n,AOA = φ_LOS,AOA + X_n · φ_n,AOA + Y_n            ... (7.5-11)
```

Ray 오프셋:
```
φ_n,m,AOA = φ_n,AOA + c_ASA · α_m                     ... (7.5-13)  (α_m: Table 7.5-3)
```

**ZOA** (Laplacian):
```
θ_n,ZOA = (2/C_θ) · (sign(U)-0.5) · sqrt(-ln(2P_n))   ... (7.5-14) (U ~ Uniform)
```

C_θ^NLOS 값: Table 7.5-4 참조. LOS인 경우:
```
C_θ = C_θ^NLOS · (1.3086 - 0.0339K - 0.0077K²)        ... (7.5-15)
```

ZOA 랜덤화:
```
θ_n,ZOA = θ_LOS,ZOA + (1-2·O2I) · X_n · θ_n,ZOA       ... (7.5-16)  (O2I=1이면 부호 반전)
θ_n,m,ZOA = θ_n,ZOA + c_ZSA · α_m                      ... (7.5-18)
```
θ_n,m,ZOA > 180°이면 360° - θ_n,m,ZOA로 래핑

**ZOD**:
```
θ_n,ZOD = θ_LOS,ZOD + X_n · θ_n,ZOD                    ... (7.5-19)
θ_n,m,ZOD = θ_n,ZOD + μ_offset,ZOD + c_ZSD · α_m       ... (7.5-20)
```

### Step 8: 클러스터 내 Ray Coupling
→ `generateSSP.cpp`

- AOD ↔ AOA 랜덤 커플링 (클러스터 n 내, 또는 sub-cluster 내)
- ZOD ↔ ZOA 동일 방식
- AOD ↔ ZOD 동일 방식

### Step 9: Cross Polarization Power Ratio (XPR)
```
κ_n,m = 10^(X_XPR / 20)                              ... (7.5-21)
```
X_XPR ~ N(μ_XPR, σ²_XPR), Table 7.5-6에서 값 참조. 각 ray·cluster 독립 생성.

**Note:** Step 1~9의 결과는 co-sited sector → 동일 UT에 대해 동일해야 함.

### Step 10: 초기 Random Phase

각 ray m, 클러스터 n에 대해 4개 편파 조합 (θθ, θφ, φθ, φφ) 각각 uniform(-π, π)에서 추출.

### Step 11: 채널 계수 생성
→ `channel.cpp`: 핵심 수식

N-2개 약한 클러스터 (n=3,...,N):
```
H_n^{us} = Σ_{m=1}^{M} [
  sqrt(P_n/M) ·
  [F_rx,u,θ(θ_n,m,ZOA, φ_n,m,AOA)]   [exp(jΦ_θθ)·sqrt(κ)   exp(jΦ_θφ)]   [F_tx,s,θ(θ_n,m,ZOD, φ_n,m,AOD)]
  [F_rx,u,φ(θ_n,m,ZOA, φ_n,m,AOA)] · [exp(jΦ_φθ)            exp(jΦ_φφ)·sqrt(κ)] · [F_tx,s,φ(θ_n,m,ZOD, φ_n,m,AOD)]
  · exp(j2π(r̂_rx · d̂_rx,u)/λ)   ← 수신 array response
  · exp(j2π(r̂_tx · d̂_tx,s)/λ)   ← 송신 array response
  · exp(j2πν_n,m·t)               ← Doppler
]                                                     ... (7.5-22)
```

2개 강한 클러스터 (n=1,2)는 3개 sub-cluster로 분할:
```
τ_{n,1} = τ_n
τ_{n,2} = τ_n + 1.28 · c_DS
τ_{n,3} = τ_n + 2.56 · c_DS                          ... (7.5-26)
```

Sub-cluster ray 매핑 (Table 7.5-5):
- k=1: rays 1~10 (power 10/20)
- k=2: rays 11~16 (power 6/20)
- k=3: rays 17~20 (power 4/20)

**LOS 채널 계수:**
```
H_LOS = [F_rx,u,θ]   [1  0]   [F_tx,s,θ]
        [F_rx,u,φ] · [0 -1] · [F_tx,s,φ]
        · exp(-j2π d_3D/λ) · array_response · Doppler   ... (7.5-29)
```

최종 채널:
```
h = sqrt(1/(1+K_R)) · h_NLOS + sqrt(K_R/(1+K_R)) · H_LOS   ... (7.5-30)
```

**Doppler:**
```
ν_n,m = (v/λ) · [sinθ_v·cosφ_v·sinθ_ZOA·cosφ_AOA + sinθ_v·sinφ_v·sinθ_ZOA·sinφ_AOA + cosθ_v·cosθ_ZOA]   ... (7.5-25)
```

### Step 12: Path Loss 및 Shadowing 적용

채널 계수에 path loss와 shadow fading 곱하기.

---

## 3. Path Loss 모델 (Section 7.4.1)

fc: 주파수 [GHz], d_2D: 수평 거리 [m], d_3D: 3D 거리 [m], h_BS/h_UT: 안테나 높이 [m]

### UMa (0.5~100 GHz)

**LOS:**
```
PL1 = 28.0 + 22·log10(d_3D) + 20·log10(fc)                     (10m ≤ d_2D ≤ d'_BP)
PL2 = 28.0 + 40·log10(d_3D) + 20·log10(fc)
      - 9·log10((d'_BP)² + (h_BS-h_UT)²)                       (d'_BP ≤ d_2D ≤ 5000m)
σ_SF = 4 dB
d'_BP = 4·h'_BS·h'_UT·fc·10⁹/c     (c = 3×10⁸ m/s)
h'_BS = h_BS - h_E,   h'_UT = h_UT - h_E
h_E = 1m (확률적: P(h_E=h) ~ uniform, h ∈ {12(n_fl-1)+1.5 : n_fl=1,...,max_fl})
```

**NLOS:**
```
PL = max(PL_UMa-LOS, PL'_UMa-NLOS)
PL'_UMa-NLOS = 13.54 + 39.08·log10(d_3D) + 20·log10(fc) - 0.6·(h_UT-1.5)
σ_SF = 6 dB

(Optional) PL = 32.4 + 20·log10(fc) + 30·log10(d_3D),  σ_SF = 7.8 dB
```

### UMi - Street Canyon (0.5~100 GHz)

**LOS:**
```
PL1 = 32.4 + 21·log10(d_3D) + 20·log10(fc)                     (10m ≤ d_2D ≤ d'_BP)
PL2 = 32.4 + 40·log10(d_3D) + 20·log10(fc)
      - 9.5·log10((d'_BP)² + (h_BS-h_UT)²)                     (d'_BP ≤ d_2D ≤ 5000m)
σ_SF = 4 dB
```

**NLOS:**
```
PL = max(PL_UMi-LOS, PL'_UMi-NLOS)
PL'_UMi-NLOS = 35.3·log10(d_3D) + 22.4 + 21.3·log10(fc) - 0.3·(h_UT-1.5)
σ_SF = 7.82 dB

(Optional) PL = 32.4 + 20·log10(fc) + 31.9·log10(d_3D),  σ_SF = 8.2 dB
```

### InH - Office (0.5~100 GHz)

**LOS:**
```
PL = 32.4 + 17.3·log10(d_3D) + 20·log10(fc)
σ_SF = 3 dB
```

**NLOS:**
```
PL = max(PL_InH-LOS, PL'_InH-NLOS)
PL'_InH-NLOS = 38.3·log10(d_3D) + 17.30 + 24.9·log10(fc)
σ_SF = 8.03 dB

(Optional) PL = 32.4 + 20·log10(fc) + 31.9·log10(d_3D),  σ_SF = 8.29 dB
```

### RMa (0.5~7 GHz)

**LOS:**
```
PL1 = 20·log10(40π·d_3D·fc/3) + min(0.03h^1.72, 10)·log10(d_3D)
      - min(0.044h^1.72, 14.77) + 0.002·log10(h)·d_3D
PL2 = PL1(d_BP) + 40·log10(d_3D/d_BP)
d_BP = 2π·h_BS·h_UT·fc·10⁹/c
σ_SF = 4 dB (d_2D < d_BP), 6 dB (d_2D ≥ d_BP)
```

**NLOS:**
```
PL = max(PL_RMa-LOS, PL'_RMa-NLOS)
PL'_RMa-NLOS = 161.04 - 7.1·log10(W) + 7.5·log10(h)
  - (24.37 - 3.7·(h/h_BS)²)·log10(h_BS)
  + (43.42 - 3.1·log10(h_BS))·(log10(d_3D) - 3)
  + 20·log10(fc) - (3.2·(log10(11.75·h_UT))² - 4.97)
σ_SF = 8 dB
```

### InF (0.5~100 GHz)

**LOS:**
```
PL = 31.84 + 21.50·log10(d_3D) + 19.00·log10(fc),  σ_SF = 4.3 dB
```

**NLOS (4 sub-scenarios):**
```
InF-SL: PL = 33.0  + 25.5·log10(d_3D) + 20·log10(fc),  σ_SF = 5.7 dB
InF-DL: PL = 18.6  + 35.7·log10(d_3D) + 20·log10(fc),  σ_SF = 7.2 dB
InF-SH: PL = 32.4  + 23.0·log10(d_3D) + 20·log10(fc),  σ_SF = 5.9 dB
InF-DH: PL = 33.63 + 21.9·log10(d_3D) + 20·log10(fc),  σ_SF = 4.0 dB
```

---

## 4. LOS 확률 모델 (Section 7.4.2)

### UMa
```
P_LOS = min(18/d_2D, 1) · (1 - exp(-d_2D/63)) + exp(-d_2D/63)     (h_UT ≤ 13m)
        1/(1 + C(d_2D, h_UT))                                       (13m < h_UT ≤ 23m)
C(d_2D, h_UT) = ((h_UT-13)/10)^1.5 · g(d_2D)
g(d_2D) = (1.25e-6)·d_2D² · exp(-d_2D/150)    (d_2D > 18m)
         = 0                                     (d_2D ≤ 18m)
```

### UMi - Street Canyon
```
P_LOS = min(18/d_2D, 1) · (1 - exp(-d_2D/36)) + exp(-d_2D/36)
```

### InH - Office
```
Mixed office:  P_LOS = { 1                       d_2D ≤ 1.2m
                        { exp(-(d_2D-1.2)/4.7)    1.2m < d_2D < 6.5m
                        { exp(-(d_2D-6.5)/32.6)·0.32    6.5m ≤ d_2D

Open office:   P_LOS = { 1                       d_2D ≤ 5m
                        { exp(-(d_2D-5)/70.8)     5m < d_2D ≤ 49m
                        { exp(-(d_2D-49)/211.7)·0.54    49m < d_2D
```

### RMa
```
P_LOS = { 1                                      d_2D ≤ 10m
         { exp(-(d_2D-10)/1000)                   10m < d_2D
```

---

## 5. O2I Penetration Loss (Section 7.4.3)

### 전체 구조
```
PL = PL_b + PL_tw + PL_in + N(0, σ_P²)
```
- PL_b = 기본 outdoor path loss
- PL_tw = 건물 외벽 관통 손실
- PL_in = 건물 내부 손실

### 재료별 관통 손실
```
L_glass    = 2 + 0.2·f    [dB]   (표준 유리, 3cm)
L_IRRglass = 25.4 + 0.11·f [dB]   (IRR 유리)
L_concrete = 5 + 4·f       [dB]   (콘크리트, 23cm)
L_plywood  = 1.03 + 0.17·f [dB]   (합판, 1.75cm)
L_wood     = 4.85 + 0.12·f [dB]   (목재, 3.3cm)
```
(f: GHz)

### O2I 모델

| 모델 | PL_tw [dB] | PL_in [dB] | σ_P [dB] |
|---|---|---|---|
| **Low-loss** | 5 - 10·log10(0.3·10^(-L_glass/10) + 0.7·10^(-L_concrete/10)) | 0.5·d_2D-in | 4.4 |
| **High-loss** | 5 - 10·log10(0.7·10^(-L_IRRglass/10) + 0.3·10^(-L_concrete/10)) | 0.5·d_2D-in | 6.5 |
| **Low-loss A** | 5 - 10·log10(0.3·10^(-L_glass/10) + 0.7·10^(-L_plywood/10)) | 0.5·d_2D-in | 4.4 |

적용:
- UMa, UMi: Low-loss + High-loss (비율은 시나리오 의존)
- Calibration (Table 7.8-1): **50% low-loss + 50% high-loss** → pathloss_model `_B` 사용

### O2I Car Penetration
```
μ = 9 dB (일반), μ = 20 dB (금속 코팅 차창)
σ_P = 5 dB
```

---

## 6. 안테나 모델링 (Section 7.3)

### BS 안테나 어레이 구조

```
(M_g, N_g, M, N, P, d_H, d_V, d_{g,H}, d_{g,V})
```
- M_g × N_g = 패널 수 (수직 × 수평)
- M × N = 패널당 엘리먼트 수 (수직 × 수평, 동일 편파)
- P = 편파 수 (1 또는 2)
- d_H, d_V = 엘리먼트 간격 [λ]
- d_{g,H}, d_{g,V} = 패널 간격 [λ]

**엘리먼트 위치 계산** (Initiallization.cpp):
```
d_tx[m][n][p][mg][ng].y = n · d_H + ng · d_{g,H}    (수평)
d_tx[m][n][p][mg][ng].z = m · d_V + mg · d_{g,V}    (수직)
```

### 단일 엘리먼트 방사 패턴 (Table 7.3-1)

**수직 컷:** (θ_3dB = 65°)
```
A_dB(θ, φ=0°) = -min{ 12·(θ/65°)², SLA_V },   SLA_V = 20 dB
```

**수평 컷:** (φ_3dB = 65°)
```
A_dB(θ=90°, φ) = -min{ 12·(φ/65°)², SLA_H },   SLA_H = 20 dB
```

**3D 패턴:**
```
A_dB(θ, φ) = -min{ -(A_dB(θ,0°) + A_dB(90°,φ)), SLA_H }
G_E,max = 8 dBi
```

### Antenna Port Mapping (Beamtilt)

```
w_m = (1/√M) · exp(-j·2π/λ · (m-1)·d_V · cos(θ_etilt))     ... (7.3-1)
```

### Polarized Antenna Modelling

**Model-2** (calibration에서 사용):
```
F_θ'(θ',φ') = sqrt(A'(θ',φ')) · cos(ζ)     ... (7.3-4)
F_φ'(θ',φ') = sqrt(A'(θ',φ')) · sin(ζ)     ... (7.3-5)
```
- ζ = 0°: 수직 편파
- ζ = ±45°: X-pol 쌍

**Model-1** (더 정확):
```
[F_θ']   [+cosψ  -sinψ]   [F_θ'']
[F_φ'] = [+sinψ  +cosψ] · [F_φ'']     ... (7.3-3)
```

### UT 안테나 (Handheld)

Directional 패턴 (Table 7.3-2):
```
θ_3dB = 125°,  SLA_V = 22.5 dB
φ_3dB = 125°,  A_max = 22.5 dB
G_E,max = 5.3 dBi
```

---

## 7. LSP 파라미터 테이블 (Table 7.5-6)

### Part 1: UMi-Street Canyon & UMa

fc: GHz, d_2D: km

#### Delay Spread (DS)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| μ_lgDS | -0.18·log10(1+fc)-7.28 | -0.22·log10(1+fc)-6.87 | -6.62 | -7.067-0.0794·log10(fc) | -6.47-0.134·log10(fc) | -6.62 |
| σ_lgDS | 0.39 | 0.19·log10(1+fc)+0.22 | 0.32 | 0.57+0.026·log10(fc) | 0.39 | 0.32 |

#### AOD Spread (ASD)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| μ_lgASD | -0.05·log10(1+fc)+1.21 | -0.24·log10(1+fc)+1.54 | 1.25 | 0.92 | 1.09 | 0.58 |
| σ_lgASD | 0.08·log10(1+fc)+0.29 | 0.10·log10(1+fc)+0.33 | 0.42 | 0.31 | 0.44 | 0.7 |

#### AOA Spread (ASA)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| μ_lgASA | -0.07·log10(1+fc)+1.66 | -0.07·log10(1+fc)+1.76 | 1.76 | 1.76 | 2.04-0.25·log10(fc) | 1.76 |
| σ_lgASA | 0.021·log10(1+fc)+0.26 | 0.05·log10(1+fc)+0.27 | 0.16 | 0.19 | 0.17-0.03·log10(fc) | 0.16 |

#### ZOA Spread (ZSA)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| μ_lgZSA | -0.11·log10(1+fc)+0.81 | -0.03·log10(1+fc)+0.92 | 1.01 | 0.96 | -0.2856·log10(fc)+1.445 | 1.01 |
| σ_lgZSA | -0.03·log10(1+fc)+0.29 | -0.05·log10(1+fc)+0.35 | 0.43 | 0.15 | 0.17 | 0.43 |

#### K-factor, XPR, Clusters

| 파라미터 | UMi LOS | UMi NLOS | UMa LOS | UMa NLOS |
|---|---|---|---|---|
| μ_K [dB] | 9 | N/A | 9 | N/A |
| σ_K [dB] | 5 | N/A | 3.5 | N/A |
| μ_XPR [dB] | 9 | 8 | 8 | 7 |
| σ_XPR [dB] | 3 | 3 | 4 | 3 |
| N_clusters | 12 | 19 | 12 | 20 |
| M_rays | 20 | 20 | 20 | 20 |
| r_τ | 3 | 2.1 | 2.5 | 2.3 |
| c_DS [ns] | 5 | 11 | max(0.25, 6.5622-3.4084·log10(fc)) | max(0.25, 6.5622-3.4084·log10(fc)) |
| c_ASD [°] | 3 | 10 | 3.58 | 1.8 |
| c_ASA [°] | 17 | 22 | 11 | 15 |
| c_ZSA [°] | 7 | 7 | 7 | 7 |
| ζ [dB] | 3 | 3 | 3 | 3 |

#### Cross-Correlations (UMi LOS / UMi NLOS / UMa LOS / UMa NLOS)

| 쌍 | UMi LOS | UMi NLOS | UMa LOS | UMa NLOS |
|---|---|---|---|---|
| ASD-DS | 0.5 | 0 | 0.4 | 0.4 |
| ASA-DS | 0.8 | 0.4 | 0.8 | 0.6 |
| ASA-SF | -0.4 | -0.4 | -0.5 | 0 |
| ASD-SF | -0.5 | 0 | -0.5 | -0.6 |
| DS-SF | -0.4 | -0.7 | -0.4 | -0.4 |
| ASD-ASA | 0.4 | 0 | 0 | 0.4 |
| ASD-K | -0.2 | N/A | 0 | N/A |
| ASA-K | -0.3 | N/A | -0.2 | N/A |
| DS-K | -0.7 | N/A | -0.4 | N/A |
| SF-K | 0.5 | N/A | 0 | N/A |
| ZSD-SF | 0 | 0 | 0 | 0 |
| ZSA-SF | 0 | 0 | -0.8 | -0.4 |
| ZSD-DS | 0 | -0.5 | -0.2 | -0.5 |
| ZSA-DS | 0.2 | 0 | 0 | 0 |
| ZSD-ASD | 0.5 | 0.5 | 0.5 | 0.5 |
| ZSA-ASD | 0.3 | 0.5 | 0 | -0.1 |
| ZSD-ASA | 0 | 0 | -0.3 | 0 |
| ZSA-ASA | 0 | 0.2 | 0.4 | 0 |
| ZSD-ZSA | 0 | 0 | 0 | 0 |

### Part 2: RMa & Indoor-Office

#### Delay Spread (DS)

| 파라미터 | RMa LOS | RMa NLOS | InH LOS | InH NLOS |
|---|---|---|---|---|
| μ_lgDS | -7.49 | -7.43 | -0.01·log10(1+fc)-7.692 | -0.28·log10(1+fc)-7.173 |
| σ_lgDS | 0.55 | 0.48 | 0.18 | 0.10·log10(1+fc)+0.055 |

#### AOD/AOA/ZOA Spread

| 파라미터 | RMa LOS | RMa NLOS | InH LOS | InH NLOS |
|---|---|---|---|---|
| μ_lgASD | 0.90 | 0.95 | 1.60 | 1.62 |
| σ_lgASD | 0.38 | 0.45 | 0.18 | 0.25 |
| μ_lgASA | 1.52 | 1.52 | -0.19·log10(1+fc)+1.781 | -0.11·log10(1+fc)+1.863 |
| σ_lgASA | 0.24 | 0.13 | 0.12·log10(1+fc)+0.119 | 0.12·log10(1+fc)+0.059 |
| μ_lgZSA | 0.47 | 0.58 | -0.26·log10(1+fc)+1.44 | -0.15·log10(1+fc)+1.387 |
| σ_lgZSA | 0.40 | 0.37 | -0.04·log10(1+fc)+0.264 | -0.09·log10(1+fc)+0.746 |

#### K-factor, XPR, Clusters

| 파라미터 | RMa LOS | RMa NLOS | InH LOS | InH NLOS |
|---|---|---|---|---|
| μ_K [dB] | 7 | N/A | 7 | N/A |
| σ_K [dB] | 4 | N/A | 4 | N/A |
| μ_XPR [dB] | 12 | 7 | 11 | 10 |
| σ_XPR [dB] | 4 | 3 | 4 | 4 |
| N_clusters | 11 | 10 | 15 | 19 |
| r_τ | 3.8 | 1.7 | 3.6 | 3 |
| c_DS [ns] | N/A | N/A | N/A | N/A |
| c_ASD [°] | 2 | 2 | 5 | 5 |
| c_ASA [°] | 3 | 3 | 8 | 11 |
| c_ZSA [°] | 3 | 3 | 9 | 9 |
| ζ [dB] | 3 | 3 | 6 | 3 |

**NOTE:** fc < 6 GHz일 때 UMa는 fc=6으로, fc < 2 GHz일 때 UMi는 fc=2로 치환

---

## 8. ZSD/ZOD 오프셋 (Table 7.5-7~12)

### UMa (Table 7.5-7)

| | LOS | NLOS |
|---|---|---|
| μ_lgZSD | max[-0.5, -2.1·(d_2D/1000) - 0.01·(h_UT-1.5) + 0.75] | max[-0.5, -2.1·(d_2D/1000) - 0.01·(h_UT-1.5) + 0.9] |
| σ_lgZSD | 0.40 | 0.49 |
| μ_offset,ZOD | 0 | e(fc) - 10^{a(fc)·log10(max(b(fc), d_2D)) + c(fc) - 0.07·(h_UT-1.5)} |

NLOS ZOD offset 계수:
```
a(fc) = 0.208·log10(fc) - 0.782
b(fc) = 25
c(fc) = -0.13·log10(fc) + 2.03
e(fc) = 7.66·log10(fc) - 5.96
```

### UMi (Table 7.5-8)

| | LOS | NLOS |
|---|---|---|
| μ_lgZSD | max[-0.21, -14.8·(d_2D/1000) + 0.01·|h_UT-h_BS| + 0.83] | max[-0.5, -3.1·(d_2D/1000) + 0.01·max(h_UT-h_BS,0) + 0.2] |
| σ_lgZSD | 0.35 | 0.35 |
| μ_offset,ZOD | 0 | -10^{-1.5·log10(max(10, d_2D)) + 3.3} |

### RMa (Table 7.5-9)

| | LOS | NLOS / O2I |
|---|---|---|
| μ_lgZSD | max[-1, -0.17·(d_2D/1000) - 0.01·(h_UT-1.5) + 0.22] | max[-1, -0.19·(d_2D/1000) - 0.01·(h_UT-1.5) + 0.28] |
| σ_lgZSD | 0.34 | 0.30 |
| μ_offset,ZOD | 0 | atan((35-3.5)/d_2D) - atan((35-1.5)/d_2D) |

### InH (Table 7.5-10)

| | LOS | NLOS |
|---|---|---|
| μ_lgZSD | -1.43·log10(1+fc) + 2.228 | 1.08 |
| σ_lgZSD | 0.13·log10(1+fc) + 0.30 | 0.36 |
| μ_offset,ZOD | 0 | 0 |

### InF (Table 7.5-11)

| | LOS | NLOS |
|---|---|---|
| μ_lgZSD | 1.35 | 1.2 |
| σ_lgZSD | 0.35 | 0.55 |
| μ_offset,ZOD | 0 | 0 |

---

## 9. Scaling Factor 테이블

### Table 7.5-2: AOA/AOD Scaling Factor (C_φ^NLOS)

| N_clusters | 4 | 5 | 6 | 7 | 8 | 10 | 11 | 12 | 14 | 15 | 16 | 19 | 20 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C_φ | 0.779 | 0.860 | 0.921 | 0.973 | 1.018 | 1.090 | 1.123 | 1.146 | 1.190 | 1.211 | 1.226 | 1.273 | 1.289 | 1.358 |

### Table 7.5-4: ZOA/ZOD Scaling Factor (C_θ^NLOS)

| N_clusters | 6 | 7 | 8 | 10 | 11 | 12 | 14 | 15 | 16 | 19 | 20 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C_θ | 0.788 | 0.847 | 0.889 | 0.957 | 1.031 | 1.104 | 1.1072 | 1.1088 | 1.1276 | 1.184 | 1.178 | 1.282 |

### Table 7.5-3: Ray Offset Angles (α_m, rms 1로 정규화)

| Ray m | α_m |
|---|---|
| 1,2 | ±0.0447 |
| 3,4 | ±0.1413 |
| 5,6 | ±0.2492 |
| 7,8 | ±0.3715 |
| 9,10 | ±0.5129 |
| 11,12 | ±0.6797 |
| 13,14 | ±0.8844 |
| 15,16 | ±1.1481 |
| 17,18 | ±1.5195 |
| 19,20 | ±2.1551 |

### Table 7.5-5: Sub-cluster 매핑 (2개 강한 클러스터)

| Sub-cluster k | Rays | Power 비율 | Delay offset |
|---|---|---|---|
| 1 | 1~10 | 10/20 | 0 |
| 2 | 11~16 | 6/20 | 1.28·c_DS |
| 3 | 17~20 | 4/20 | 2.56·c_DS |

---

## 10. Calibration 파라미터 (Section 7.8)

### Table 7.8-1: Large Scale Calibration

| Parameter | Values |
|---|---|
| Scenarios | UMa, UMi-Street Canyon, Indoor-office (open office) |
| Sectorization | 3 sectors: 30°, 150°, 270° |
| BS antenna | M_g=N_g=1, (M,N,P)=(10,1,1), d_V=0.5λ |
| BS port mapping | 10 elements → 1 CRS port |
| BS downtilt | 102° (UMa, UMi), 110° (Indoor) |
| Antenna virtualization | DFT precoding (TR 36.897) |
| BS Tx power | 44 dBm (UMi@6G), 49 dBm (UMa@6G), 35 dBm (@30/70G), 24 dBm (Indoor) |
| Bandwidth | 20 MHz (6G), 100 MHz (30/70G) |
| UT antenna | 1 element (vertical), Isotropic |
| Handover margin | 0 dB |
| UT distribution | TR36.873 (UMa/UMi 3D), uniform (Indoor, min 0m) |
| UT attachment | Pathloss-based (LOS angle) |
| UT noise figure | 9 dB |
| Fast fading | Not modelled |
| O2I penetration | 50% low-loss + 50% high-loss |
| Carrier frequency | 6, 30, 70 GHz |
| Wrapping | Geographical distance based (mandatory) |
| Metrics | 1) Coupling loss CDF, 2) Geometry CDF (with/without noise) |

### Table 7.8-2: Full Calibration

| Parameter | Values |
|---|---|
| Scenarios | UMa, UMi-Street Canyon, Indoor-office (open office) |
| Carrier frequency | 6, 30, 60, 70 GHz |
| Bandwidth | 20 MHz (6G), 100 MHz (30/60/70G) |
| BS Tx power | Same as Table 7.8-1 |
| **BS Config 1** | M=4, N=4, P=2, M_g=1, N_g=2, d_H=d_V=0.5λ, d_{H,g}=d_{V,g}=2.5λ → metrics 1,2,3 |
| BS Config 1 port | 16 elem/pol/panel → 1 CRS port, panning (0,0)° |
| **BS Config 2** | M_g=N_g=1, M=N=2, P=1 → metrics 1,2,4 |
| BS Config 2 port | Each element → 1 CRS port |
| UT antenna | M_g=N_g=1, M=N=1, P=2 |
| UT attachment | RSRP (formula) from CRS port 0 |
| Polarized antenna | Model-2 (TR36.873) |
| UT orientation | Ω_α ~ U[0,360°], Ω_β = 90°, Ω_γ = 0° |
| UT antenna pattern | Isotropic |
| Polarization | P=2 → BS: X-pol (±45°), UT: X-pol (0/+90°) |
| Metrics | 1) Coupling loss, 2) Wideband SIR, 3) DS/AS CDFs, 4) PRB SVD CDFs |
