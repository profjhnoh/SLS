# Handheld UT 안테나 모델 구현 상세 문서

> **기준 스펙**: 3GPP TR 38.901 V17, Section 7.8, Table 7.8-2A Config B
> **구현 일자**: 2026-02-19
> **대상 시나리오**: Dense Urban 7 GHz Calibration

---

## 1. 개요

### 1.1 배경

기존 시뮬레이터의 UE 안테나 모델은 두 가지:
- **Omni** (`ue_antenna_element_gain = 0`): 등방성 안테나, 편파만 구분
- **Directional** (`ue_antenna_element_gain > 0`): 단일 boresight 방향의 방향성 안테나

TR 38.901 Section 7.8에서는 **handheld UT**를 정의한다. 이 모델에서 UE는 $15\text{cm} \times 7\text{cm}$ 직사각형 디바이스이며, 8개의 후보 안테나 위치가 있고, 각 안테나는 **고유한 boresight 방향**과 **편파 방향**을 가진다.

### 1.2 Config B 사양 (7 GHz Calibration)

| 항목 | 값 |
|------|-----|
| 안테나 포트 수 | 4 |
| 안테나 위치 (1-based) | 1, 7, 3, 5 |
| 편파 | Single polarization ($\zeta = 0$) |
| Element pattern | Table 7.3-2 ($\theta_{3\text{dB}} = \varphi_{3\text{dB}} = 125°$, $G_{\max} = 5.3$ dBi) |
| UE orientation | $\alpha \sim U[0°, 360°]$, $\beta = 45°$, $\gamma = 0°$ |
| RSRP 계산 | 4포트 평균 |

### 1.3 핵심 차이: 2단계 좌표 회전

기존 방향성 안테나는 **1단계 회전**만 수행한다:

$$
\text{기존}: \quad F'(\theta', \varphi') \xrightarrow{\text{Stage 2: } (\Omega_\alpha, \Omega_\beta, \Omega_\gamma)} F(\theta, \varphi) \text{ in GCS}
$$

Handheld 모델은 안테나마다 고유한 방향이 있으므로 **2단계 회전**이 필요하다:

$$
F'(\theta'', \varphi'') \xrightarrow{\text{Stage 1: } (\alpha_u, 90°, 0°)} F'(\theta', \varphi') \text{ in UT LCS} \xrightarrow{\text{Stage 2: } (\Omega_\alpha, \Omega_\beta, \Omega_\gamma)} F(\theta, \varphi) \text{ in GCS}
$$

여기서:
- **Stage 1**: 안테나 개별 좌표계 → UT 본체 좌표계 (안테나별 상수)
- **Stage 2**: UT 본체 좌표계 → 글로벌 좌표계 (UE 방향에 따라 달라짐)

---

## 2. 디바이스 물리 모델

### 2.1 안테나 위치

$15\text{cm} \times 7\text{cm}$ 직사각형의 중심을 원점으로, 8개 후보 위치는 꼭짓점과 변의 중점이다:

```
     7 ─────── 6 ─────── 5
     │                     │
     │                     │
     8       (center)      4
     │                     │
     │                     │
     1 ─────── 2 ─────── 3

     ← 15 cm (x축) →
     ↕ 7 cm (y축)
```

| 위치 | 좌표 $(x, y, z)$ [m] | 설명 |
|------|----------------------|------|
| 1 | $(-0.075, -0.035, 0)$ | 좌하 꼭짓점 |
| 2 | $(0, -0.035, 0)$ | 하변 중점 |
| 3 | $(+0.075, -0.035, 0)$ | 우하 꼭짓점 |
| 4 | $(+0.075, 0, 0)$ | 우변 중점 |
| 5 | $(+0.075, +0.035, 0)$ | 우상 꼭짓점 |
| 6 | $(0, +0.035, 0)$ | 상변 중점 |
| 7 | $(-0.075, +0.035, 0)$ | 좌상 꼭짓점 |
| 8 | $(-0.075, 0, 0)$ | 좌변 중점 |

### 2.2 안테나별 Euler 각도 ($\alpha_u$)

각 안테나의 boresight는 디바이스 중심에서 안테나 위치를 향한다. 이 방향의 방위각이 $\alpha_u$이다:

$$
\alpha_u = \text{atan2}(x_{\text{pos}}, -y_{\text{pos}})
$$

| 위치 | $\alpha_u$ (rad) | 근사값 (°) |
|------|------------------|-----------|
| 1 | $\text{atan2}(-0.075, 0.035)$ | $\approx -65°$ |
| 2 | $\text{atan2}(0, 0.035)$ | $0°$ |
| 3 | $\text{atan2}(0.075, 0.035)$ | $\approx +65°$ |
| 4 | $\text{atan2}(0.075, 0)$ | $90°$ |
| 5 | $\text{atan2}(0.075, -0.035)$ | $\approx +115°$ |
| 6 | $\text{atan2}(0, -0.035)$ | $180°$ |
| 7 | $\text{atan2}(-0.075, -0.035)$ | $\approx -115°$ |
| 8 | $\text{atan2}(-0.075, 0)$ | $-90°$ |

모든 안테나에 대해 $\beta_u = 90°$, $\gamma_u = 0°$ (상수).

**Config B에서 사용하는 포트**: 위치 1, 7, 3, 5 → 디바이스 네 꼭짓점

---

## 3. 수학적 모델: 2단계 좌표 회전

### 3.1 좌표 변환 함수 (기존 코드 재사용)

기존 코드의 `Get_LCS_theta()`와 `Get_LCS_pi()`는 Euler 각도 $(\alpha, \beta, \gamma)$를 이용해 GCS 방향을 LCS 방향으로 변환한다:

$$
\theta' = \arccos\!\big[\cos\beta\cos\gamma\cos\theta + (\sin\beta\cos\gamma\cos(\varphi-\alpha) - \sin\gamma\sin(\varphi-\alpha))\sin\theta\big]
$$

$$
\varphi' = \arg\!\big[(\cos\beta\sin\theta\cos(\varphi-\alpha) - \sin\beta\cos\theta) + j(\cos\beta\sin\gamma\cos\theta + (\sin\beta\sin\gamma\cos(\varphi-\alpha) + \cos\gamma\sin(\varphi-\alpha))\sin\theta)\big]
$$

### 3.2 2단계 역변환 절차

GCS 방향 $(\theta, \varphi)$가 주어졌을 때:

**Step 1 — Stage 2 역변환** (GCS → UT LCS):
$$
\theta' = \texttt{Get\_LCS\_theta}(\Omega_\alpha, \Omega_\beta, \Omega_\gamma, \theta, \varphi)
$$
$$
\varphi' = \texttt{Get\_LCS\_pi}(\Omega_\alpha, \Omega_\beta, \Omega_\gamma, \theta, \varphi)
$$

**Step 2 — Stage 1 역변환** (UT LCS → 안테나 기준 좌표계):
$$
\theta'' = \texttt{Get\_LCS\_theta}(\alpha_u, \beta_u, \gamma_u, \theta', \varphi')
$$
$$
\varphi'' = \texttt{Get\_LCS\_pi}(\alpha_u, \beta_u, \gamma_u, \theta', \varphi')
$$

### 3.3 Element Pattern (Table 7.3-2)

안테나 기준 좌표계에서의 방향 $(\theta'', \varphi'')$를 이용하여:

$$
A_V = -\min\!\left(12\left(\frac{\theta'' - 90°}{125°}\right)^2, \; 22.5\right) \quad \text{[dB]}
$$

$$
A_H = -\min\!\left(12\left(\frac{\varphi''}{125°}\right)^2, \; 22.5\right) \quad \text{[dB]}
$$

$$
A(\theta'', \varphi'') = G_{\max} - \min\!\big(-(A_V + A_H), \; 22.5\big) \quad \text{[dBi]}
$$

여기서 $G_{\max} = 5.3$ dBi, $\theta_{3\text{dB}} = \varphi_{3\text{dB}} = 125°$, $A_{\max} = 22.5$ dB.

진폭으로 변환:
$$
A_{\text{lin}} = \sqrt{10^{A/10}}
$$

### 3.4 필드 패턴 및 2단계 순변환

Single polarization ($\zeta = 0$)이므로 안테나 기준 좌표계에서:
$$
F'_\theta = A_{\text{lin}} \cdot \cos(\zeta) = A_{\text{lin}}, \quad F'_\varphi = A_{\text{lin}} \cdot \sin(\zeta) = 0
$$

**Stage 1 순변환** (안테나 기준 → UT LCS):
$$
\texttt{LCS\_to\_GCS}(\alpha_u, \beta_u, \gamma_u, \; \theta', \varphi', \; F'_\theta, F'_\varphi) \;\rightarrow\; (F'_{\theta,\text{UT}}, \; F'_{\varphi,\text{UT}})
$$

**Stage 2 순변환** (UT LCS → GCS):
$$
\texttt{LCS\_to\_GCS}(\Omega_\alpha, \Omega_\beta, \Omega_\gamma, \; \theta, \varphi, \; F'_{\theta,\text{UT}}, F'_{\varphi,\text{UT}}) \;\rightarrow\; (F_\theta, \; F_\varphi) \text{ in GCS}
$$

최종 $(F_\theta, F_\varphi)$는 이 안테나 포트의 GCS field pattern이다.

---

## 4. RSRP 계산

### 4.1 포트별 RSRP

각 포트 $u$에 대한 RSRP:

$$
\text{RSRP}_u = \underbrace{\frac{K_R}{K_R+1} \left|AF_{\text{tx,LOS}}\right|^2 \left|\mathbf{F}_{\text{rx},u}^T \cdot \mathbf{H}_{\text{LOS}} \cdot \mathbf{F}_{\text{tx}}\right|^2}_{\text{LOS component}} + \underbrace{\sum_{n=1}^{N} \frac{P_n}{M_n(K_R+1)} \sum_{m=1}^{M_n} \left|AF_{\text{tx}}\right|^2 \left|\mathbf{F}_{\text{rx},u}^T \cdot \mathbf{H}_{n,m} \cdot \mathbf{F}_{\text{tx}}\right|^2}_{\text{NLOS component}}
$$

여기서:
- $AF_{\text{tx}}$: BS TX array factor (BS 빔 가중치에 의존)
- $\mathbf{F}_{\text{rx},u}$: 포트 $u$의 GCS field pattern (2단계 회전 결과)
- $\mathbf{H}_{n,m}$: 편파 행렬 (XPR, 랜덤 위상 포함)

**핵심**: 기존 코드에서 $\mathbf{F}_{\text{rx}}$는 모든 포트에 동일했지만, handheld에서는 **포트마다 다르다**.

### 4.2 포트 평균 RSRP

$$
\text{RSRP} = \frac{1}{N_{\text{ports}}} \sum_{u=1}^{N_{\text{ports}}} \text{RSRP}_u
$$

Config B에서 $N_{\text{ports}} = 4$.

---

## 5. 코드 변경 사항

### 5.1 전역변수 선언

**`h/common.h`** — extern 선언 추가:
```cpp
extern int handheld_mode;            // 0=off, 1=handheld
extern int handheld_num_ports;       // 활성 포트 수 (e.g., 4)
extern int handheld_port_indices[8]; // 1-based 안테나 위치 인덱스
extern Real handheld_beta_deg;       // UT β 각도 (도)
```

**`h/Initiallization.h`** — 실제 정의:
```cpp
int handheld_mode = 0;
int handheld_num_ports = 4;
int handheld_port_indices[8] = {1,7,3,5,0,0,0,0};
Real handheld_beta_deg = 45.0;
```

### 5.2 CFG 파라미터 파싱

**`src/setSimulParam.cpp`** — `Set_simul_param()` 함수 내:

```cpp
handheld_mode       = int(Get_parameter(infile, "handheld_mode", 0));
handheld_beta_deg   = Real(Get_parameter(infile, "handheld_beta", 45.0));
handheld_num_ports  = int(Get_parameter(infile, "handheld_num_ports", 4));
handheld_port_indices[0] = int(Get_parameter(infile, "handheld_port_1", 1));
handheld_port_indices[1] = int(Get_parameter(infile, "handheld_port_2", 7));
handheld_port_indices[2] = int(Get_parameter(infile, "handheld_port_3", 3));
handheld_port_indices[3] = int(Get_parameter(infile, "handheld_port_4", 5));
```

CFG 파일 예시:
```
handheld_mode         1
handheld_beta         45
handheld_num_ports    4
handheld_port_1       1
handheld_port_2       7
handheld_port_3       3
handheld_port_4       5

MS_num_of_ant_elements_with_same_pol_each_col(M)  4
MS_num_of_ant_elements_columns(N)                  1
MS_antenna_polarization(P)                         1
MS_Mp                                              4
MS_Np                                              1
ue_antenna_element_gain                            5.3
```

### 5.3 UE 방향 초기화

**`src/Initiallization.cpp`** — `Set_antenna_location_vector()` 함수:

**변경 전**: `ue_antenna_element_gain == 0`이면 $\alpha=0$, 아니면 $\alpha \sim U[0, 2\pi]$, $\beta=\gamma=0$.

**변경 후**: handheld 모드일 때 별도 분기:
```cpp
if (handheld_mode) {
    ms[ms_idx].alpha = 360 * randnum.u() * (pi / 180.);  // U[0, 2π]
    ms[ms_idx].beta  = handheld_beta_deg * (pi / 180.);   // 45° → 0.785 rad
    ms[ms_idx].gamma = 0.;
}
```

기존 방향성 안테나: $\beta = 0$ (수직 자세)
Handheld: $\beta = 45°$ (기울어진 자세, 3GPP 규격)

### 5.4 안테나 위치 (`d_rx`) 초기화

**변경 전**: 정규 그리드 (`m * dV`, `n * dH`)로 안테나 원소 배치 후 UE 방향으로 회전.

**변경 후**: handheld 모드일 때 물리적 디바이스 좌표를 직접 사용:

```cpp
const LOCATION3D handheld_positions[8] = {
    {-0.075, -0.035, 0},  // 위치 1
    { 0,     -0.035, 0},  // 위치 2
    { 0.075, -0.035, 0},  // 위치 3
    { 0.075,  0,     0},  // 위치 4
    { 0.075,  0.035, 0},  // 위치 5
    { 0,      0.035, 0},  // 위치 6
    {-0.075,  0.035, 0},  // 위치 7
    {-0.075,  0,     0},  // 위치 8
};

for (int port = 0; port < handheld_num_ports; port++) {
    int ant_idx = handheld_port_indices[port] - 1;
    ms[ue_idx].d_rx[port][0][0][0][0] = Transpose_LCS_to_GCS_location(
        ms[ue_idx].alpha, ms[ue_idx].beta, ms[ue_idx].gamma,
        handheld_positions[ant_idx]);
}
```

`d_rx[port][0][0][0][0]` 매핑:
- 인덱스 `[M][N][P][Mg][Ng]`에서 `M = port index` (0~3)
- `MS_M = 4`, `MS_N = 1`, `MS_P = 1` → 4개 포트, 각각 단일 원소

### 5.5 핵심 변경: `Get_UE_antenna_pattern()` — 2단계 회전

**`src/Link.cpp`** — 함수 시그니처 변경:

```cpp
// 변경 전
Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS,
    int ms_idx, int sector_index,
    Real &F_theta_GCS_P1, Real &F_pi_GCS_P1,
    Real &F_theta_GCS_P2, Real &F_pi_GCS_P2);

// 변경 후 (default parameter 추가)
Real Get_UE_antenna_pattern(int P, Real theta_GCS, Real pi_GCS,
    int ms_idx, int sector_index,
    Real &F_theta_GCS_P1, Real &F_pi_GCS_P1,
    Real &F_theta_GCS_P2, Real &F_pi_GCS_P2,
    int port_idx = -1);  // ← 새 파라미터
```

`port_idx = -1`: 기존 동작 (handheld 아님)
`port_idx >= 0`: handheld 포트 인덱스 (0-based)

**handheld 분기 구현** (함수 최상단에 추가):

```cpp
if (handheld_mode && port_idx >= 0)
{
    // (1) 상수 테이블: 안테나별 α_u
    static const Real HANDHELD_ALPHA_U[8] = {
        atan2(-0.075,  0.035),  // pos 1: ≈ -65°
        atan2( 0.0,    0.035),  // pos 2:    0°
        atan2( 0.075,  0.035),  // pos 3: ≈ +65°
        atan2( 0.075,  0.0  ),  // pos 4:   90°
        atan2( 0.075, -0.035),  // pos 5: ≈ +115°
        atan2( 0.0,   -0.035),  // pos 6:  180°
        atan2(-0.075, -0.035),  // pos 7: ≈ -115°
        atan2(-0.075,  0.0  ),  // pos 8:  -90°
    };

    int ant_idx = handheld_port_indices[port_idx] - 1;
    Real alpha_u = HANDHELD_ALPHA_U[ant_idx];
    Real beta_u  = π / 2;   // 90°
    Real gamma_u = 0;

    // (2) Stage 2 역변환: GCS → UT LCS
    Real θ' = Get_LCS_theta(Ω_α, Ω_β, Ω_γ, θ, φ);
    Real φ' = Get_LCS_pi   (Ω_α, Ω_β, Ω_γ, θ, φ);

    // (3) Stage 1 역변환: UT LCS → 안테나 기준 좌표계
    Real θ'' = Get_LCS_theta(α_u, β_u, γ_u, θ', φ');
    Real φ'' = Get_LCS_pi   (α_u, β_u, γ_u, θ', φ');

    // (4) Element pattern 계산 (Table 7.3-2)
    // → combined_gain [dBi], A = sqrt(10^(gain/10))

    // (5) Stage 1 순변환: F'(θ'', φ'') → F'(θ', φ') in UT LCS
    LCS_Antenna_field_to_GCS_antenna_pattern(α_u, β_u, γ_u,
        θ', φ', F'_θ, F'_φ, F_θ_UT, F_φ_UT);

    // (6) Stage 2 순변환: F(θ', φ') → F(θ, φ) in GCS
    LCS_Antenna_field_to_GCS_antenna_pattern(Ω_α, Ω_β, Ω_γ,
        θ, φ, F_θ_UT, F_φ_UT, F_θ_GCS, F_φ_GCS);

    // (7) Panel 2 = 0 (single pol, single panel)
    return combined_gain;
}
```

**기존 TYPE==11/12/13 분기는 변경 없이 유지.**

### 5.6 호출부 조건 변경: omni 판별 기준

Handheld 모드에서는 `ue_antenna_element_gain = 5.3` (≠ 0)으로 설정하므로, 기존 분기 `if (ue_antenna_element_gain == 0)`는 자연스럽게 BF 경로를 탄다.

추가로, 아래 3개 함수에서 omni 조건을 보강하여 안전하게 처리:

| 함수 | 파일 | 변경 |
|------|------|------|
| `Get_MS_antgain()` | Link.cpp:2490 | `ue_antenna_element_gain == 0` → `ue_antenna_element_gain == 0 && !handheld_mode` |
| `Get_RX_SmallScale_antgain()` | Link.cpp:1993, 2065 | 동일 패턴 (InH, DU/Rural 각각) |

이 함수들에서 directional 경로를 탈 때 `port_idx`를 전달:
```cpp
Get_UE_antenna_pattern(P, v_angle_theta, h_angle_pi, ms_idx, sector_idx,
    F_theta_GCS_P1, F_pi_GCS_P1, F_theta_GCS_P2, F_pi_GCS_P2,
    handheld_mode ? M : -1);  // M = port index (루프 변수)
```

### 5.7 `find_best_tx_beam()` — 포트별 ray_power

**문제**: 기존 코드는 Phase 1에서 `ray_power[polRx][n][m]`을 한 번만 계산한다. Handheld에서는 포트마다 element pattern이 다르므로 포트별로 따로 계산해야 한다.

**해결**: 포트별 배열 `hh_ray_power[port][n][m]` 추가.

```
Phase 1 (beam-independent, 1회 계산):
┌─────────────────────────────────────────────────┐
│  for each ray (n, m):                           │
│    BS TX pattern: F_tx(θ, φ) — 공통             │
│    Random phases, XPR — 공통                     │
│                                                  │
│    if (handheld):                                │
│      for each port u = 0..3:                     │
│        F_rx,u = Get_UE_antenna_pattern(port=u)   │
│        hh_ray_power[u][n][m] = |F_rx,u·H·F_tx|² │
│    else:                                         │
│        ... (기존 코드)                            │
└─────────────────────────────────────────────────┘

Phase 2 (beam search, TX 빔마다 반복):
┌─────────────────────────────────────────────────┐
│  for each TX beam (a, z):                        │
│    for each ray (n, m):                          │
│      AF_tx² = |Σ w·exp(j·phase)|²               │
│      Σ_port hh_ray_power[port][n][m]             │
│      alpha += P_n/M_n × AF_tx² × Σ_port         │
│                                                  │
│    RSRP = alpha / num_ports                      │
└─────────────────────────────────────────────────┘
```

LOS 성분도 동일하게 포트별 처리: `hh_los_ray_power[port]`.

### 5.8 `Get_RSRP()` — 변경 불필요

`ue_antenna_element_gain = 5.3` (≠ 0)이므로 기존 else 분기 (UE beamforming 경로)를 자동으로 탄다.

이 경로에서:
- `MS_Mp = 4` → `m` 루프가 0~3 (= port index)
- `Get_MS_antgain(m, ...)` → `Get_UE_antenna_pattern(..., port_idx=m)`
- `Get_RX_SmallScale_antgain(m, ...)` → 동일
- 최종 `alpha / (M*N*P) = alpha / (4*1*1)` = 포트 평균

기존 BF 경로와 정확히 일치하므로 추가 수정 불필요.

---

## 6. 데이터 흐름 요약

```
┌──────────────────────────────────────────────────────────────────┐
│                        CFG 파일 설정                              │
│  handheld_mode = 1, handheld_beta = 45                           │
│  MS_M = 4, MS_N = 1, MS_P = 1, MS_Mp = 4, MS_Np = 1            │
│  ue_antenna_element_gain = 5.3                                   │
│  handheld_port_1..4 = 1, 7, 3, 5                                │
└──────────────────────┬───────────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────────┐
│               UE 초기화 (Set_antenna_location_vector)             │
│                                                                   │
│  α ~ U[0, 2π],  β = 45° (0.785 rad),  γ = 0                    │
│                                                                   │
│  d_rx[0][0][0][0][0] = Rotate(α,β,γ, pos_1)  ← 위치 1           │
│  d_rx[1][0][0][0][0] = Rotate(α,β,γ, pos_7)  ← 위치 7           │
│  d_rx[2][0][0][0][0] = Rotate(α,β,γ, pos_3)  ← 위치 3           │
│  d_rx[3][0][0][0][0] = Rotate(α,β,γ, pos_5)  ← 위치 5           │
└──────────────────────┬───────────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────────┐
│              RSRP 계산 (find_best_tx_beam)                        │
│                                                                   │
│  Phase 1: for each ray (n,m), for each port u:                   │
│    GCS → UT LCS → 안테나 기준 좌표계 (2단계 역변환)                │
│    Element pattern A(θ'', φ'') 계산                               │
│    안테나 기준 → UT LCS → GCS (2단계 순변환)                       │
│    → hh_ray_power[u][n][m]                                       │
│                                                                   │
│  Phase 2: for each TX beam:                                      │
│    RSRP = (1/4) × Σ_port Σ_n Σ_m P_n/M_n × |AF_tx|²            │
│            × hh_ray_power[port][n][m]                             │
└──────────────────────────────────────────────────────────────────┘
```

---

## 7. 수정 파일 목록

| 파일 | 변경 내용 |
|------|----------|
| `h/common.h` | `handheld_mode`, `handheld_num_ports`, `handheld_port_indices`, `handheld_beta_deg` extern 선언 |
| `h/Initiallization.h` | 위 변수들의 실제 정의 (초기값 포함) |
| `src/setSimulParam.cpp` | CFG 파라미터 파싱 + 로그 출력 |
| `src/Initiallization.cpp` | UE orientation ($\beta=45°$) + d_rx 물리적 위치 설정 |
| `src/Link.cpp` | `Get_UE_antenna_pattern()`: port_idx 파라미터 + 2단계 회전 구현 |
| `src/Link.cpp` | `Get_MS_antgain()`: handheld omni 분기 조건 보강 |
| `src/Link.cpp` | `Get_RX_SmallScale_antgain()`: 동일 (2곳) |
| `src/Link.cpp` | `find_best_tx_beam()`: 포트별 hh_ray_power 계산 + handheld Phase 2 |
| `src/channel.cpp` | forward declaration에 port_idx 추가 |

---

## 8. 기존 코드와의 호환성

- `handheld_mode = 0` (기본값)일 때 **모든 변경 사항이 비활성화**됨
- `port_idx` 파라미터는 default value `-1`이므로 기존 호출부 수정 불필요
- 기존 omni/directional 경로는 `handheld_mode = 0`일 때 100% 동일하게 동작
- 전역변수 초기값이 `handheld_mode = 0`이므로 cfg에 해당 파라미터가 없어도 안전

---

## 9. 검증 계획

1. **Regression**: `handheld_mode = 0` → 기존 결과와 bit-exact 동일
2. **포트별 패턴 확인**: 4개 포트의 $(F_\theta, F_\varphi)$가 서로 다른지 로그 출력
3. **UE α 회전 확인**: $\alpha$가 달라질 때 편파 패턴이 올바르게 회전하는지 확인
4. **Coupling loss CDF**: Table 7.8-2A Config B 기준 곡선과 비교
