# 3GPP TR 38.901 Quick Reference for Tbps_DLSLS Simulator

> **Source:** 3GPP TR 38.901 V19.1.0 (2025-09) — Study on channel model for frequencies from 0.5 to 100 GHz
> **Purpose:** 시뮬레이터 개발 시 빠른 참조용. 전체 원문은 `doc/38901-j10.docx` 참고.

---

## 목차

1. [코드-Spec 매핑](#1-코드-spec-매핑)
2. [채널 생성 절차 (Section 7.5)](#2-채널-생성-절차-section-75) — Step 1~12 상세
3. [Path Loss 모델 (Section 7.4.1)](#3-path-loss-모델-section-741) — UMa, UMi, InH, RMa, InF, SMa
4. [LOS 확률 모델 (Section 7.4.2)](#4-los-확률-모델-section-742)
5. [O2I Penetration Loss (Section 7.4.3)](#5-o2i-penetration-loss-section-743)
6. [안테나 모델링 (Section 7.3)](#6-안테나-모델링-section-73)
7. [LSP 파라미터 테이블 (Table 7.5-6)](#7-lsp-파라미터-테이블-table-756) — Part 1~4: UMi/UMa/RMa/InH/InF/SMa + Cross-Correlations + Correlation Distance
8. [ZSD/ZOD 오프셋 (Table 7.5-7~12)](#8-zsdzod-오프셋-table-757~12)
9. [Scaling Factor 테이블](#9-scaling-factor-테이블) — Table 7.5-2,3,4,5
10. [Calibration 파라미터 (Section 7.8)](#10-calibration-파라미터-section-78)
11. [추가 모델링 컴포넌트 (Section 7.6)](#11-추가-모델링-컴포넌트-section-76) — 산소 흡수, 공간 일관성, Blockage
12. [RSRP 계산 공식 (TR 36.873)](#12-rsrp-계산-공식-tr-36873)
13. [Annex A: Circular Angle Spread 계산](#13-annex-a-circular-angle-spread-계산)

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

> 채널 계수 생성은 Figure 7.5-1의 절차를 따르며, 하향링크(downlink)를 가정한다. 상향링크의 경우 arrival/departure 파라미터를 교환한다.
> GCS에서 zenith각 $\theta=0°$는 천정, $\theta=90°$는 수평선이다. 구면 단위벡터 $\hat{\theta}$, $\hat{\varphi}$는 전파 방향 기준으로 정의된다.

### Table 7.5-1: GCS 표기법 (Notations)

| Parameter | Notation |
|---|---|
| LOS AOD | $\phi_{\text{LOS,AOD}}$ |
| LOS AOA | $\phi_{\text{LOS,AOA}}$ |
| LOS ZOD | $\theta_{\text{LOS,ZOD}}$ |
| LOS ZOA | $\theta_{\text{LOS,ZOA}}$ |
| Cluster $n$ AOA | $\phi_{n,\text{AOA}}$ |
| Cluster $n$ AOD | $\phi_{n,\text{AOD}}$ |
| Cluster $n$, ray $m$ AOA | $\phi_{n,m,\text{AOA}}$ |
| Cluster $n$, ray $m$ AOD | $\phi_{n,m,\text{AOD}}$ |
| Cluster $n$ ZOA | $\theta_{n,\text{ZOA}}$ |
| Cluster $n$ ZOD | $\theta_{n,\text{ZOD}}$ |
| Cluster $n$, ray $m$ ZOA | $\theta_{n,m,\text{ZOA}}$ |
| Cluster $n$, ray $m$ ZOD | $\theta_{n,m,\text{ZOD}}$ |
| Rx antenna element $u$, $\theta$ field pattern | $F_{rx,u,\theta}$ |
| Rx antenna element $u$, $\varphi$ field pattern | $F_{rx,u,\varphi}$ |
| Tx antenna element $s$, $\theta$ field pattern | $F_{tx,s,\theta}$ |
| Tx antenna element $s$, $\varphi$ field pattern | $F_{tx,s,\varphi}$ |

---

### Step 1: 환경/레이아웃/안테나 설정
→ `Initiallization.cpp`: `InitializeSystem()`, `Initialdrop()`

**(a)** 시나리오 선택: UMa, UMi-Street Canyon, RMa, InH-Office, InF 중 하나를 선택한다. GCS(Global Coordinate System)를 정의하고 zenith각 $\theta$, azimuth각 $\varphi$, 구면 단위벡터 $\hat{\theta}$, $\hat{\varphi}$를 정의한다.
- **Note:** RMa는 7 GHz까지, 나머지는 100 GHz까지 적용 가능.

**(b)** BS 및 UT 개수를 설정한다.

**(c)** BS와 UT의 3D 위치를 결정하고, GCS에서 각 BS-UT 쌍의 LOS 방향을 계산한다:
- LOS AOD ($\phi_{\text{LOS,AOD}}$), LOS ZOD ($\theta_{\text{LOS,ZOD}}$)
- LOS AOA ($\phi_{\text{LOS,AOA}}$), LOS ZOA ($\theta_{\text{LOS,ZOA}}$)

**(d)** BS 및 UT 안테나 field pattern $F_{rx}$, $F_{tx}$를 GCS 기준으로 설정한다. Array geometry를 정의한다.

**(e)** BS 및 UT orientation을 GCS에 대해 설정한다:
- BS: $\Omega_{BS,\alpha}$ (bearing), $\Omega_{BS,\beta}$ (downtilt), $\Omega_{BS,\gamma}$ (slant)
- UT: $\Omega_{UT,\alpha}$ (bearing), $\Omega_{UT,\beta}$ (downtilt), $\Omega_{UT,\gamma}$ (slant)

**(f)** UT의 속도와 이동 방향을 GCS에서 설정한다.

**(g)** 시스템 중심 주파수 $f_c$와 대역폭 $B$를 설정한다.

> **Note:** Wrapping 사용 시, BS/site의 각 wrapping copy는 별도의 BS/site로 취급하여 채널을 생성한다.

---

### Step 2: 전파 조건 할당 (Large Scale Parameters 시작)
→ `channel.cpp`: LOS/NLOS 판정

Table 7.4.2-1에 따라 각 BS-UT 링크에 전파 조건(LOS/NLOS)을 할당한다.
- **서로 다른 BS-UT 링크의 전파 조건은 uncorrelated이다.**
- Indoor/outdoor 상태를 각 UT에 할당한다. **한 UT의 모든 링크는 동일한 indoor/outdoor 상태를 가진다.**

---

### Step 3: Path Loss 계산
→ `Link.cpp`: `Get_CouplingLoss()`

Table 7.4.1-1의 수식으로 모델링 대상인 각 BS-UT 링크의 path loss를 계산한다.

> **Note:** $PL = PL_b + PL_{tw} + PL_{in} + \mathcal{N}(0, \sigma_P^2)$ (O2I의 경우), $PL_b$에서 $d_{3D}$는 $d_{3D\text{-out}}$으로 대체.

---

### Step 4: Large Scale Parameter (LSP) 생성
→ `generateLSP.cpp`: `Generate_LSP()`

Table 7.5-6의 파라미터를 사용하여 다음 7개 LSP를 생성한다:
- **DS** (Delay Spread), **ASD** (AOD Spread), **ASA** (AOA Spread)
- **ZSD** (ZOD Spread), **ZSA** (ZOA Spread)
- **K** (Ricean K-factor, LOS에서만), **SF** (Shadow Fading)

**생성 절차:**
1. Table 7.5-6의 cross-correlation 행렬을 **Cholesky 분해**하여 square root matrix를 구한다.
2. 7개 독립 가우시안 랜덤 벡터에 Cholesky matrix를 곱하여 상관된 LSP 벡터를 생성한다.
3. LSP 벡터 순서: $\mathbf{s}_M = [s_{SF},\; s_K,\; s_{DS},\; s_{ASD},\; s_{ASA},\; s_{ZSD},\; s_{ZSA}]^T$
4. 각 LSP에 대해 log-normal 변환: $X = 10^{(\mu_{\lg X} + \sigma_{\lg X} \cdot s_X)}$ (단, K와 SF는 dB 도메인에서 직접)

**규칙:**
- 서로 다른 BS-UT 링크의 LSP는 uncorrelated
- **Co-sited sector → 동일 UT에 대해 동일 LSP** (같은 사이트의 다른 섹터)
- 서로 다른 층(floor)의 UT 간 LSP는 uncorrelated

**제한:**

$$
\text{ASA} = \min(\text{ASA},\; 104°), \quad
\text{ASD} = \min(\text{ASD},\; 104°), \quad
\text{ZSA} = \min(\text{ZSA},\; 52°), \quad
\text{ZSD} = \min(\text{ZSD},\; 52°)
$$

**주파수 치환 규칙:**
- UMa: $f_c < 6$ GHz이면 $f_c = 6$으로 치환하여 LSP 결정
- UMi: $f_c < 2$ GHz이면 $f_c = 2$로 치환하여 LSP 결정
- InH: $f_c < 6$ GHz이면 $f_c = 6$으로 치환하여 LSP 결정

### Step 5: 클러스터 Delay 생성 (Small Scale Parameters 시작)
→ `generateSSP.cpp`

Table 7.5-6의 delay distribution에 따라 delay를 추출한다. 지수 분포(exponential delay distribution)의 경우:

$$
\tau'_n = -r_\tau \cdot DS \cdot \ln(X_n), \quad X_n \sim \text{Uniform}(0,1), \quad n = 1,\ldots,N \tag{7.5-1}
$$

여기서 $r_\tau$는 delay distribution proportionality factor (Table 7.5-6), $N$은 클러스터 수이다.

최소 delay를 빼고 오름차순 정렬하여 정규화:

$$
\tau_n = \text{sort}\!\left(\tau'_n - \min(\tau')\right) \tag{7.5-2}
$$

**LOS인 경우:** LOS peak 추가에 의한 delay spread 변화를 보상하기 위한 추가 스케일링이 필요하다.
K-factor 의존 스케일링 상수:

$$
C_\tau = 0.7705 - 0.0433 K + 0.0002 K^2 + 0.000017 K^3 \tag{7.5-3}
$$

여기서 $K$ [dB]는 Step 4에서 생성한 Ricean K-factor이다.

스케일된 delay:

$$
\tau_n^{\text{LOS}} = \tau_n / C_\tau \tag{7.5-4}
$$

> **중요:** 스케일된 delay $\tau_n^{\text{LOS}}$는 **클러스터 power 생성(Step 6)에는 사용하지 않는다!**
> Step 6에서는 스케일 전의 $\tau_n$ (Eq. 7.5-2)를 사용한다.

---

### Step 6: 클러스터 Power 생성

단일 기울기 지수형 PDP(Power Delay Profile)를 가정하여 클러스터 파워를 계산한다.
지수 분포 delay distribution의 경우:

$$
P'_n = \exp\!\left(-\tau_n \cdot \frac{r_\tau - 1}{r_\tau \cdot DS}\right) \cdot 10^{-\zeta_n / 10} \tag{7.5-5}
$$

여기서 $\zeta_n$은 per-cluster shadowing [dB] (Table 7.5-6의 $\zeta$ 값, 표준 분포 $\mathcal{N}(0, \zeta^2)$).

정규화하여 전체 클러스터 파워 합 = 1:

$$
P_n = \frac{P'_n}{\sum_{i} P'_i} \tag{7.5-6}
$$

**LOS인 경우:** 첫 번째 클러스터에 specular component 추가:

$$
P_{1,\text{LOS}} = \frac{K_R}{K_R + 1} \tag{7.5-7}
$$

$$
P_n = \frac{1}{K_R + 1} \cdot P_n + \delta(n-1) \cdot \frac{K_R}{K_R + 1} \tag{7.5-8}
$$

여기서 $\delta(\cdot)$는 Dirac delta 함수, $K_R$은 Step 4의 K-factor (linear scale).
**이 power 값은 Eq. (7.5-9)와 (7.5-14)에서만 사용하고, Eq. (7.5-22)에서는 사용하지 않는다.**

각 클러스터 내 ray의 power: $P_n / M$ ($M$ = rays per cluster, Table 7.5-6).

**클러스터 제거:** Eq. (7.5-6) 기준 최대 클러스터 파워 대비 **-25 dB** 미만인 클러스터를 제거한다. 제거 후 scaling factor는 변경하지 않는다.

---

### Step 7: 도래각/출발각 생성 (Azimuth + Elevation)

#### 7-1. AOA 생성 (Wrapped Gaussian)

모든 클러스터의 azimuth PAS(Power Angular Spectrum)을 wrapped Gaussian으로 모델링한다.
Inverse Gaussian 함수로 AOA를 결정:

$$
\phi'_{n,\text{AOA}} = \frac{2}{C_\phi} \cdot \frac{\text{ASA}}{1.4} \cdot \sqrt{-\ln\!\left(\frac{P_n}{\max(P_n)}\right)} \tag{7.5-9}
$$

$C_\phi$ 정의:

$$
C_\phi = C_\phi^{\text{NLOS}} \cdot \left(1.1035 - 0.028 K - 0.002 K^2 + 0.0001 K^3\right) \quad \text{(LOS)} \tag{7.5-10}
$$

$$
C_\phi = C_\phi^{\text{NLOS}} \quad \text{(NLOS)}
$$

여기서 $C_\phi^{\text{NLOS}}$는 클러스터 수에 따른 scaling factor (**Table 7.5-2** 참조), $K$는 [dB] 단위.

**Table 7.5-2: AOA/AOD Scaling Factor ($C_\phi^{\text{NLOS}}$)**

| $N_{\text{clusters}}$ | 4 | 5 | 6 | 7 | 8 | 10 | 11 | 12 | 14 | 15 | 16 | 19 | 20 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $C_\phi$ | 0.779 | 0.860 | 0.921 | 0.973 | 1.018 | 1.090 | 1.123 | 1.146 | 1.190 | 1.211 | 1.226 | 1.273 | 1.289 | 1.358 |

부호 할당 및 랜덤 오프셋 추가:

$$
\phi_{n,\text{AOA}} = X_n \cdot \phi'_{n,\text{AOA}} + Y_n + \phi_{\text{LOS,AOA}} \tag{7.5-11}
$$

여기서 $X_n \in \{+1, -1\}$ uniform, $Y_n \sim \mathcal{N}(0, (\text{ASA}/7)^2)$.

**LOS인 경우:** 첫 번째 클러스터를 LOS 방향으로 강제:

$$
\phi_{n,\text{AOA}} = \phi_{n,\text{AOA}} - \left(X_1 \cdot \phi'_{1,\text{AOA}} + Y_1 - \phi_{\text{LOS,AOA}}\right) \tag{7.5-12}
$$

Per-ray 오프셋 추가 (**Table 7.5-3** 참조):

$$
\phi_{n,m,\text{AOA}} = \phi_{n,\text{AOA}} + c_{\text{ASA}} \cdot \alpha_m \tag{7.5-13}
$$

여기서 $c_{\text{ASA}}$는 cluster ASA (Table 7.5-6), $\alpha_m$은 Table 7.5-3의 ray offset angle.

**Table 7.5-3: Ray Offset Angles ($\alpha_m$, rms 1로 정규화)**

| Ray $m$ | $\alpha_m$ |
|---|---|
| 1, 2 | $\pm 0.0447$ |
| 3, 4 | $\pm 0.1413$ |
| 5, 6 | $\pm 0.2492$ |
| 7, 8 | $\pm 0.3715$ |
| 9, 10 | $\pm 0.5129$ |
| 11, 12 | $\pm 0.6797$ |
| 13, 14 | $\pm 0.8844$ |
| 15, 16 | $\pm 1.1481$ |
| 17, 18 | $\pm 1.5195$ |
| 19, 20 | $\pm 2.1551$ |

#### 7-2. AOD 생성

AOA와 동일한 절차를 따르되, ASA → ASD, $c_{\text{ASA}}$ → $c_{\text{ASD}}$, $\phi_{\text{LOS,AOA}}$ → $\phi_{\text{LOS,AOD}}$로 교체한다.

#### 7-3. ZOA 생성 (Laplacian)

모든 클러스터의 zenith PAS를 **Laplacian** 분포로 모델링한다.
Inverse Laplacian 함수로 ZOA를 결정:

$$
\theta'_{n,\text{ZOA}} = -\frac{\text{ZSA}}{C_\theta} \cdot \ln\!\left(\frac{P_n}{\max(P_n)}\right) \tag{7.5-14}
$$

$C_\theta$ 정의:

$$
C_\theta = C_\theta^{\text{NLOS}} \cdot \left(1.3086 + 0.0339 K - 0.0077 K^2 + 0.0002 K^3\right) \quad \text{(LOS)} \tag{7.5-15}
$$

$$
C_\theta = C_\theta^{\text{NLOS}} \quad \text{(NLOS)}
$$

여기서 $C_\theta^{\text{NLOS}}$는 **Table 7.5-4** 참조.

**Table 7.5-4: ZOA/ZOD Scaling Factor ($C_\theta^{\text{NLOS}}$)**

| $N_{\text{clusters}}$ | 6 | 7 | 8 | 10 | 11 | 12 | 14 | 15 | 16 | 19 | 20 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $C_\theta$ | 0.788 | 0.847 | 0.889 | 0.957 | 1.031 | 1.104 | 1.1072 | 1.1088 | 1.1276 | 1.184 | 1.178 | 1.282 |

부호 할당 및 랜덤 오프셋 추가:

$$
\theta_{n,\text{ZOA}} = X_n \cdot \theta'_{n,\text{ZOA}} + Y_n + \theta_{\text{LOS,ZOA}} \quad \text{(non-O2I)} \tag{7.5-16}
$$

$$
\theta_{n,\text{ZOA}} = X_n \cdot \theta'_{n,\text{ZOA}} + Y_n + 90° \quad \text{(O2I, 부호 반전)}
$$

여기서 $X_n \in \{+1, -1\}$ uniform, $Y_n \sim \mathcal{N}(0, (\text{ZSA}/7)^2)$.
- **O2I 링크:** $(1-2 \cdot \text{O2I})$ 인자 적용 → O2I=1이면 부호가 반전됨.

**LOS인 경우:** Eq. (7.5-16)을 Eq. (7.5-17)로 대체하여 첫 클러스터를 LOS 방향으로 강제:

$$
\theta_{n,\text{ZOA}} = \theta_{n,\text{ZOA}} - \left(X_1 \cdot \theta'_{1,\text{ZOA}} + Y_1 - \theta_{\text{LOS,ZOA}}\right) \tag{7.5-17}
$$

Per-ray 오프셋:

$$
\theta_{n,m,\text{ZOA}} = \theta_{n,\text{ZOA}} + c_{\text{ZSA}} \cdot \alpha_m \tag{7.5-18}
$$

$c_{\text{ZSA}}$는 cluster ZSA (Table 7.5-6).
> **래핑:** $\theta_{n,m,\text{ZOA}} \in [0°, 360°]$로 가정. **$\theta_{n,m,\text{ZOA}} > 180°$이면 $360° - \theta_{n,m,\text{ZOA}}$로 설정.**

#### 7-4. ZOD 생성

ZOA와 동일한 절차를 따르되, Eq. (7.5-16) 대신:

$$
\theta_{n,\text{ZOD}} = X_n \cdot \theta'_{n,\text{ZOD}} + \theta_{\text{LOS,ZOD}} \tag{7.5-19}
$$

($Y_n$ 랜덤 오프셋 없음, O2I 부호 반전 없음)

Per-ray 오프셋에 $\mu_{\text{offset,ZOD}}$ 추가 (Tables 7.5-7~12):

$$
\theta_{n,m,\text{ZOD}} = \theta_{n,\text{ZOD}} + \mu_{\text{offset,ZOD}} + \frac{3}{8} \cdot 10^{\mu_{\lg\text{ZSD}}} \cdot \alpha_m \tag{7.5-20}
$$

> **Note:** ZOD의 per-ray offset에는 $c_{\text{ZSD}}$가 아니라 $\frac{3}{8} \cdot 10^{\mu_{\lg\text{ZSD}}}$를 사용한다.
> $\mu_{\lg\text{ZSD}}$는 Table 7.5-7~12의 ZSD log-normal 분포의 평균값이다.

**LOS인 경우:** ZOA와 동일하게 Eq. (7.5-17) 절차를 따른다.

---

### Step 8: 클러스터 내 Ray Coupling
→ `generateSSP.cpp`

클러스터 $n$ 내에서 (또는 2개 강한 클러스터의 sub-cluster 내에서, Step 11 참조) 다음을 랜덤하게 커플링한다:
1. **AOD ↔ AOA**: $\phi_{n,m,\text{AOD}}$를 $\phi_{n,m,\text{AOA}}$에 랜덤 매핑
2. **ZOD ↔ ZOA**: $\theta_{n,m,\text{ZOD}}$를 $\theta_{n,m,\text{ZOA}}$에 동일한 방식으로 랜덤 매핑
3. **AOD ↔ ZOD**: $\phi_{n,m,\text{AOD}}$를 $\theta_{n,m,\text{ZOD}}$에 동일한 방식으로 랜덤 매핑

---

### Step 9: Cross Polarization Power Ratio (XPR)

각 클러스터 $n$의 각 ray $m$에 대해 XPR을 log-Normal 분포에서 독립적으로 생성한다:

$$
\kappa_{n,m} = 10^{X_{n,m} / 10} \tag{7.5-21}
$$

여기서 $X_{n,m} \sim \mathcal{N}(\mu_{\text{XPR}}, \sigma^2_{\text{XPR}})$, $\mu_{\text{XPR}}$과 $\sigma_{\text{XPR}}$은 Table 7.5-6에서 참조한다.

> **Note:** $X_{n,m}$은 각 ray와 각 cluster에 대해 **독립적으로** 추출한다.

> **★ Step 1~9의 결과는 co-sited sector → 동일 UT에 대해 동일해야 한다.**

### Step 10: 초기 Random Phase 생성 (Coefficient Generation 시작)

각 ray $m$, 클러스터 $n$에 대해 4개 편파 조합의 초기 위상을 추출한다:

$$
\Phi^{\theta\theta}_{n,m},\; \Phi^{\theta\varphi}_{n,m},\; \Phi^{\varphi\theta}_{n,m},\; \Phi^{\varphi\varphi}_{n,m} \sim \text{Uniform}(-\pi, \pi)
$$

각 위상은 독립적으로 생성된다.

---

### Step 11: 채널 계수 생성
→ `channel.cpp`, `channel_update.cpp`: 핵심 수식

#### 11-1. 약한 클러스터 ($N-2$개, $n = 3, 4, \ldots, N$)

각 수신 안테나 $u$, 송신 안테나 $s$ 쌍에 대한 NLOS 채널 계수:

$$
H_n^{\text{NLOS},us}(t) = \sqrt{\frac{P_n}{M}} \sum_{m=1}^{M}
\begin{bmatrix}
F_{rx,u,\theta}(\theta_{n,m,\text{ZOA}}, \phi_{n,m,\text{AOA}}) \\
F_{rx,u,\varphi}(\theta_{n,m,\text{ZOA}}, \phi_{n,m,\text{AOA}})
\end{bmatrix}^T
\begin{bmatrix}
e^{j\Phi^{\theta\theta}_{n,m}} & \sqrt{\kappa_{n,m}^{-1}}\, e^{j\Phi^{\theta\varphi}_{n,m}} \\
\sqrt{\kappa_{n,m}^{-1}}\, e^{j\Phi^{\varphi\theta}_{n,m}} & e^{j\Phi^{\varphi\varphi}_{n,m}}
\end{bmatrix}
\begin{bmatrix}
F_{tx,s,\theta}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}}) \\
F_{tx,s,\varphi}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}})
\end{bmatrix}
\cdot e^{j \frac{2\pi}{\lambda_0} (\hat{r}_{rx,n,m} \cdot \bar{d}_{rx,u})}
\cdot e^{j \frac{2\pi}{\lambda_0} (\hat{r}_{tx,n,m} \cdot \bar{d}_{tx,s})}
\cdot e^{j 2\pi \nu_{n,m} t}
\tag{7.5-22}
$$

여기서:
- $F_{rx,u,\theta}$, $F_{rx,u,\varphi}$: 수신 안테나 엘리먼트 $u$의 GCS field pattern (Clause 7.1 변환 포함)
- $F_{tx,s,\theta}$, $F_{tx,s,\varphi}$: 송신 안테나 엘리먼트 $s$의 GCS field pattern
- $\kappa_{n,m}$: cross polarisation power ratio (linear scale, Step 9)
- $\lambda_0$: 캐리어 주파수 파장
- **편파 무시 시:** $2\times2$ 편파 행렬을 스칼라 $e^{j\Phi^{\theta\theta}_{n,m}}$로 대체하고 수직 편파 field pattern만 적용

**구면 단위벡터:**

$$
\hat{r}_{rx,n,m} = \begin{bmatrix}
\sin\theta_{n,m,\text{ZOA}} \cos\phi_{n,m,\text{AOA}} \\
\sin\theta_{n,m,\text{ZOA}} \sin\phi_{n,m,\text{AOA}} \\
\cos\theta_{n,m,\text{ZOA}}
\end{bmatrix} \tag{7.5-23}
$$

$$
\hat{r}_{tx,n,m} = \begin{bmatrix}
\sin\theta_{n,m,\text{ZOD}} \cos\phi_{n,m,\text{AOD}} \\
\sin\theta_{n,m,\text{ZOD}} \sin\phi_{n,m,\text{AOD}} \\
\cos\theta_{n,m,\text{ZOD}}
\end{bmatrix} \tag{7.5-24}
$$

**Doppler 주파수:**

$$
\nu_{n,m} = \frac{1}{\lambda_0} \hat{r}_{rx,n,m}^T \cdot \bar{v}
= \frac{v}{\lambda_0} \left[
\sin\theta_v \cos\phi_v \sin\theta_{n,m,\text{ZOA}} \cos\phi_{n,m,\text{AOA}}
+ \sin\theta_v \sin\phi_v \sin\theta_{n,m,\text{ZOA}} \sin\phi_{n,m,\text{AOA}}
+ \cos\theta_v \cos\theta_{n,m,\text{ZOA}}
\right] \tag{7.5-25}
$$

여기서 $\bar{v} = (v, \theta_v, \phi_v)$는 UT 속도 벡터.

#### 11-2. 강한 클러스터 (2개, $n = 1, 2$) → Sub-cluster 분할

2개 가장 강한 클러스터($n=1,2$)의 20개 ray를 3개 sub-cluster로 분할한다.

**Sub-cluster delay:**

$$
\tau_{n,1} = \tau_n, \quad
\tau_{n,2} = \tau_n + 1.28 \cdot c_{DS}, \quad
\tau_{n,3} = \tau_n + 2.56 \cdot c_{DS} \tag{7.5-26}
$$

여기서 $c_{DS}$는 Table 7.5-6의 cluster delay spread.
> **$c_{DS}$가 N/A인 경우:** 3.91 ns를 사용한다 (legacy 동작: sub-cluster delay 5 ns, 10 ns).

**Table 7.5-5: Sub-cluster 매핑**

| Sub-cluster $k$ | Ray 번호 | Power 비율 | Delay offset |
|---|---|---|---|
| 1 | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 | 10/20 | 0 |
| 2 | 11, 12, 13, 14, 15, 16 | 6/20 | $1.28 \cdot c_{DS}$ |
| 3 | 17, 18, 19, 20 | 4/20 | $2.56 \cdot c_{DS}$ |

**채널 임펄스 응답:**

$$
H_n^{us}(\tau, t) = \sum_{i=1}^{3} H_{n,i}^{\text{NLOS},us}(t) \cdot \delta(\tau - \tau_{n,i}) \tag{7.5-27}
$$

여기서 $H_{n,i}^{\text{NLOS},us}$는 Eq. (7.5-22)와 동일하되, sub-cluster $i$에 속한 ray만 합산:

$$
H_{n,i}^{\text{NLOS},us}(t) = \text{Eq.(7.5-22) with } m \in \text{sub-cluster } i \tag{7.5-28}
$$

#### 11-3. LOS 채널 계수

LOS 조건에서 직접 경로의 채널 계수:

$$
H_1^{\text{LOS},us}(t) =
\begin{bmatrix}
F_{rx,u,\theta}(\theta_{\text{LOS,ZOA}}, \phi_{\text{LOS,AOA}}) \\
F_{rx,u,\varphi}(\theta_{\text{LOS,ZOA}}, \phi_{\text{LOS,AOA}})
\end{bmatrix}^T
\begin{bmatrix}
1 & 0 \\
0 & -1
\end{bmatrix}
\begin{bmatrix}
F_{tx,s,\theta}(\theta_{\text{LOS,ZOD}}, \phi_{\text{LOS,AOD}}) \\
F_{tx,s,\varphi}(\theta_{\text{LOS,ZOD}}, \phi_{\text{LOS,AOD}})
\end{bmatrix}
\cdot e^{-j \frac{2\pi}{\lambda_0} d_{3D}}
\cdot e^{j \frac{2\pi}{\lambda_0} (\hat{r}_{\text{LOS,rx}} \cdot \bar{d}_{rx,u})}
\cdot e^{j \frac{2\pi}{\lambda_0} (\hat{r}_{\text{LOS,tx}} \cdot \bar{d}_{tx,s})}
\cdot e^{j 2\pi \nu_{\text{LOS}} t}
\tag{7.5-29}
$$

> **Note:** LOS 경로의 편파 행렬은 $\begin{bmatrix}1 & 0 \\ 0 & -1\end{bmatrix}$이다 (cross-pol 성분 없음).

#### 11-4. 최종 채널 임펄스 응답 (LOS)

NLOS 채널에 LOS 성분을 추가하고, K-factor에 따라 스케일링:

$$
H^{us}(\tau, t) = \sqrt{\frac{1}{K_R+1}} \cdot H^{\text{NLOS},us}(\tau, t) + \sqrt{\frac{K_R}{K_R+1}} \cdot H_1^{\text{LOS},us}(t) \cdot \delta(\tau - \tau_1) \tag{7.5-30}
$$

여기서 $K_R$은 Step 4에서 생성한 Ricean K-factor (linear scale).

---

### Step 12: Path Loss 및 Shadowing 적용

최종 채널 계수에 path loss와 shadow fading을 곱한다:

$$
H_{\text{final}} = H \cdot 10^{-(PL + SF)/20}
$$

여기서 $PL$은 Step 3의 path loss [dB], $SF$는 Step 4의 shadow fading [dB].
> **Note:** SF의 부호 정의: 양(+)의 SF는 path loss 모델 예측보다 UT에서 더 많은 수신 전력을 의미한다.

---

## 3. Path Loss 모델 (Section 7.4.1)

$f_c$: 주파수 [GHz], $d_{2D}$: 수평 거리 [m], $d_{3D}$: 3D 거리 [m], $h_{BS}$/$h_{UT}$: 안테나 높이 [m]

### UMa (0.5~100 GHz)

**LOS:**

$$
PL_1 = 28.0 + 22 \log_{10}(d_{3D}) + 20 \log_{10}(f_c) \quad (10\text{m} \le d_{2D} \le d'_{BP})
$$

$$
PL_2 = 28.0 + 40 \log_{10}(d_{3D}) + 20 \log_{10}(f_c) - 9 \log_{10}\!\left((d'_{BP})^2 + (h_{BS}-h_{UT})^2\right) \quad (d'_{BP} \le d_{2D} \le 5000\text{m})
$$

$$
\sigma_{SF} = 4 \text{ dB}
$$

$$
d'_{BP} = \frac{4 \, h'_{BS} \, h'_{UT} \, f_c \times 10^9}{c}, \quad c = 3 \times 10^8 \text{ m/s}
$$

$$
h'_{BS} = h_{BS} - h_E, \quad h'_{UT} = h_{UT} - h_E
$$

$h_E = 1$ m (확률적: $P(h_E = h)$ ~ uniform, $h \in \{12(n_{fl}-1)+1.5 : n_{fl}=1,\ldots,\text{max\_fl}\}$)

**NLOS:**

$$
PL = \max(PL_{\text{UMa-LOS}},\; PL'_{\text{UMa-NLOS}})
$$

$$
PL'_{\text{UMa-NLOS}} = 13.54 + 39.08 \log_{10}(d_{3D}) + 20 \log_{10}(f_c) - 0.6(h_{UT}-1.5)
$$

$$
\sigma_{SF} = 6 \text{ dB}
$$

(Optional) $PL = 32.4 + 20\log_{10}(f_c) + 30\log_{10}(d_{3D})$, $\sigma_{SF} = 7.8$ dB

### UMi - Street Canyon (0.5~100 GHz)

**LOS:**

$$
PL_1 = 32.4 + 21 \log_{10}(d_{3D}) + 20 \log_{10}(f_c) \quad (10\text{m} \le d_{2D} \le d'_{BP})
$$

$$
PL_2 = 32.4 + 40 \log_{10}(d_{3D}) + 20 \log_{10}(f_c) - 9.5 \log_{10}\!\left((d'_{BP})^2 + (h_{BS}-h_{UT})^2\right) \quad (d'_{BP} \le d_{2D} \le 5000\text{m})
$$

$$
\sigma_{SF} = 4 \text{ dB}
$$

**NLOS:**

$$
PL = \max(PL_{\text{UMi-LOS}},\; PL'_{\text{UMi-NLOS}})
$$

$$
PL'_{\text{UMi-NLOS}} = 22.4 + 35.3 \log_{10}(d_{3D}) + 21.3 \log_{10}(f_c) - 0.3(h_{UT}-1.5)
$$

$$
\sigma_{SF} = 7.82 \text{ dB}
$$

(Optional) $PL = 32.4 + 20\log_{10}(f_c) + 31.9\log_{10}(d_{3D})$, $\sigma_{SF} = 8.2$ dB

### InH - Office (0.5~100 GHz)

**LOS:**

$$
PL = 32.4 + 17.3 \log_{10}(d_{3D}) + 20 \log_{10}(f_c), \quad \sigma_{SF} = 3 \text{ dB}
$$

**NLOS:**

$$
PL = \max(PL_{\text{InH-LOS}},\; PL'_{\text{InH-NLOS}})
$$

$$
PL'_{\text{InH-NLOS}} = 17.30 + 38.3 \log_{10}(d_{3D}) + 24.9 \log_{10}(f_c), \quad \sigma_{SF} = 8.03 \text{ dB}
$$

(Optional) $PL = 32.4 + 20\log_{10}(f_c) + 31.9\log_{10}(d_{3D})$, $\sigma_{SF} = 8.29$ dB

### RMa (0.5~7 GHz)

**LOS:**

$$
PL_1 = 20\log_{10}\!\left(\frac{40\pi \, d_{3D} \, f_c}{3}\right) + \min(0.03 h^{1.72}, 10)\log_{10}(d_{3D}) - \min(0.044 h^{1.72}, 14.77) + 0.002 \log_{10}(h) \cdot d_{3D}
$$

$$
PL_2 = PL_1(d_{BP}) + 40\log_{10}(d_{3D}/d_{BP})
$$

$$
d_{BP} = \frac{2\pi \, h_{BS} \, h_{UT} \, f_c \times 10^9}{c}
$$

$$
\sigma_{SF} = 4 \text{ dB } (d_{2D} < d_{BP}), \quad 6 \text{ dB } (d_{2D} \ge d_{BP})
$$

**NLOS:**

$$
PL = \max(PL_{\text{RMa-LOS}},\; PL'_{\text{RMa-NLOS}})
$$

$$
PL'_{\text{RMa-NLOS}} = 161.04 - 7.1\log_{10}(W) + 7.5\log_{10}(h) - (24.37 - 3.7(h/h_{BS})^2)\log_{10}(h_{BS})
$$
$$
+ (43.42 - 3.1\log_{10}(h_{BS}))(\log_{10}(d_{3D}) - 3) + 20\log_{10}(f_c) - (3.2(\log_{10}(11.75 h_{UT}))^2 - 4.97)
$$

$$
\sigma_{SF} = 8 \text{ dB}
$$

### InF (0.5~100 GHz)

**LOS:**

$$
PL = 31.84 + 21.50 \log_{10}(d_{3D}) + 19.00 \log_{10}(f_c), \quad \sigma_{SF} = 4.3 \text{ dB}
$$

**NLOS (4 sub-scenarios):**

$$
\text{InF-SL}: \; PL = 33.0 + 25.5 \log_{10}(d_{3D}) + 20 \log_{10}(f_c), \quad \sigma_{SF} = 5.7 \text{ dB}
$$

$$
\text{InF-DL}: \; PL = 18.6 + 35.7 \log_{10}(d_{3D}) + 20 \log_{10}(f_c), \quad \sigma_{SF} = 7.2 \text{ dB}
$$

$$
\text{InF-SH}: \; PL = 32.4 + 23.0 \log_{10}(d_{3D}) + 20 \log_{10}(f_c), \quad \sigma_{SF} = 5.9 \text{ dB}
$$

$$
\text{InF-DH}: \; PL = 33.63 + 21.9 \log_{10}(d_{3D}) + 20 \log_{10}(f_c), \quad \sigma_{SF} = 4.0 \text{ dB}
$$

### SMa (0.5~37 GHz, new in Rel-19)

**LOS:**

$$
PL_1 = 26.0 + 22 \log_{10}(d_{3D}) + 20 \log_{10}(f_c) \quad (10\text{m} \le d_{2D} \le d'_{BP})
$$

$$
PL_2 = PL_1(d'_{BP}) + 40\log_{10}(d_{3D}/d'_{BP}) \quad (d'_{BP} \le d_{2D} \le 5000\text{m})
$$

$$
\sigma_{SF} = 4 \text{ dB}
$$

$$
d'_{BP} = \frac{4 \, h'_{BS} \, h'_{UT} \, f_c \times 10^9}{c}
$$

**NLOS:**

$$
PL = \max(PL_{\text{SMa-LOS}},\; PL'_{\text{SMa-NLOS}})
$$

$$
PL'_{\text{SMa-NLOS}} = 13.54 + 39.08 \log_{10}(d_{3D}) + 20 \log_{10}(f_c) - 0.6(h_{UT}-1.5)
$$

$$
\sigma_{SF} = 6 \text{ dB}
$$

> **Note:** $f_c$ 정규화 = 1 GHz, 거리 정규화 = 1 m. $0.5 < f_c < 37$ GHz.

---

## 4. LOS 확률 모델 (Section 7.4.2)

### UMa

$$
P_{\text{LOS}} =
\begin{cases}
\min\!\left(\dfrac{18}{d_{2D}}, 1\right) \cdot \left(1 - e^{-d_{2D}/63}\right) + e^{-d_{2D}/63} & h_{UT} \le 13\text{m} \\[6pt]
\dfrac{1}{1 + C(d_{2D}, h_{UT})} & 13\text{m} < h_{UT} \le 23\text{m}
\end{cases}
$$

$$
C(d_{2D}, h_{UT}) = \left(\frac{h_{UT}-13}{10}\right)^{1.5} \cdot g(d_{2D})
$$

$$
g(d_{2D}) =
\begin{cases}
1.25 \times 10^{-6} \cdot d_{2D}^2 \cdot e^{-d_{2D}/150} & d_{2D} > 18\text{m} \\
0 & d_{2D} \le 18\text{m}
\end{cases}
$$

### UMi - Street Canyon

$$
P_{\text{LOS}} = \min\!\left(\frac{18}{d_{2D}}, 1\right) \cdot \left(1 - e^{-d_{2D}/36}\right) + e^{-d_{2D}/36}
$$

### InH - Office

$$
\text{Mixed office:} \quad P_{\text{LOS}} =
\begin{cases}
1 & d_{2D} \le 1.2\text{m} \\
e^{-(d_{2D}-1.2)/4.7} & 1.2\text{m} < d_{2D} < 6.5\text{m} \\
0.32 \cdot e^{-(d_{2D}-6.5)/32.6} & 6.5\text{m} \le d_{2D}
\end{cases}
$$

$$
\text{Open office:} \quad P_{\text{LOS}} =
\begin{cases}
1 & d_{2D} \le 5\text{m} \\
e^{-(d_{2D}-5)/70.8} & 5\text{m} < d_{2D} \le 49\text{m} \\
0.54 \cdot e^{-(d_{2D}-49)/211.7} & 49\text{m} < d_{2D}
\end{cases}
$$

### RMa

$$
P_{\text{LOS}} =
\begin{cases}
1 & d_{2D} \le 10\text{m} \\
e^{-(d_{2D}-10)/1000} & 10\text{m} < d_{2D}
\end{cases}
$$

### InF (Indoor Factory)

$$
\text{InF-SL:} \quad P_{\text{LOS}} = e^{-d_{2D} / k_{SL}}
$$

$$
\text{InF-DL:} \quad P_{\text{LOS}} = e^{-d_{2D} / k_{DL}}
$$

$$
\text{InF-SH:} \quad P_{\text{LOS}} = 1 - \left(1 - e^{-d_{2D}/k_{SH1}}\right)\left(1 + \frac{d_{2D}}{k_{SH2}} \cdot e^{-d_{2D}/k_{SH2}}\right)
$$

$$
\text{InF-DH:} \quad P_{\text{LOS}} = 1 - \left(1 - e^{-d_{2D}/k_{DH1}}\right)\left(1 + \frac{d_{2D}}{k_{DH2}} \cdot e^{-d_{2D}/k_{DH2}}\right)
$$

$$
\text{InF-HH:} \quad P_{\text{LOS}} = 1 \quad \text{(항상 LOS)}
$$

> Note: $k$ 값들은 clutter 밀도($r$), clutter 높이($h_c$), BS/UT 높이에 의존한다. 상세 수식은 Table 7.4.2-1 참조.

### SMa (Suburban Macro)

$$
P_{\text{LOS}} = e^{-d_{2D} / 63} \quad (d_{\text{out}})
$$

> Note: outdoor 부분 거리에 적용. 상세 수식은 Table 7.4.2-1의 SMa 항 참조.

---

## 5. O2I Penetration Loss (Section 7.4.3)

### 전체 구조

$$
PL = PL_b + PL_{tw} + PL_{in} + \mathcal{N}(0, \sigma_P^2) \tag{7.4-2}
$$

- $PL_b$ = 기본 outdoor path loss
- $PL_{tw}$ = 건물 외벽 관통 손실
- $PL_{in}$ = 건물 내부 손실

### 재료별 관통 손실

$$
L_{\text{glass}} = 2 + 0.2 f \;\text{[dB]} \quad \text{(표준 유리, 3cm)}
$$

$$
L_{\text{IRRglass}} = 23 + 0.3 f \;\text{[dB]} \quad \text{(IRR 유리)}
$$

$$
L_{\text{concrete}} = 5 + 4 f \;\text{[dB]} \quad \text{(콘크리트, 23cm)}
$$

$$
L_{\text{plywood}} = 1.03 + 0.17 f \;\text{[dB]} \quad \text{(합판, 1.75cm)}
$$

$$
L_{\text{wood}} = 4.85 + 0.12 f \;\text{[dB]} \quad \text{(목재, 3.3cm)}
$$

($f$: GHz)

### O2I 모델

| 모델 | $PL_{tw}$ [dB] | $PL_{in}$ [dB] | $\sigma_P$ [dB] |
|---|---|---|---|
| **Low-loss** | $5 - 10\log_{10}(0.3 \cdot 10^{-L_{\text{glass}}/10} + 0.7 \cdot 10^{-L_{\text{concrete}}/10})$ | $0.5 \cdot d_{\text{2D-in}}$ | 4.4 |
| **High-loss** | $5 - 10\log_{10}(0.7 \cdot 10^{-L_{\text{IRRglass}}/10} + 0.3 \cdot 10^{-L_{\text{concrete}}/10})$ | $0.5 \cdot d_{\text{2D-in}}$ | 6.5 |
| **Low-loss A** | $5 - 10\log_{10}(0.3 \cdot 10^{-L_{\text{glass}}/10} + 0.7 \cdot 10^{-L_{\text{plywood}}/10})$ | $0.5 \cdot d_{\text{2D-in}}$ | 4.4 |

$d_{\text{2D,in}}$ 생성:
- UMa, UMi-Street Canyon: $d_{\text{2D,in}} = \min(25 U_1, 25 U_2)$, $U_1, U_2 \sim \text{Uniform}(0,1)$ 독립
- RMa: $d_{\text{2D,in}} = \min(10 U_1, 10 U_2)$, $U_1, U_2 \sim \text{Uniform}(0,1)$ 독립
- UT별로 독립 생성

적용:
- UMa, UMi: Low-loss + High-loss (비율은 시나리오 의존)
- Calibration (Table 7.8-1): **50% low-loss + 50% high-loss** → pathloss_model `_B` 사용

### O2I Car Penetration

$$
PL = PL_b + \mathcal{N}(\mu, \sigma_P^2) \tag{7.4-4}
$$

$\mu = 9$ dB (일반), $\mu = 20$ dB (금속 코팅 차창), $\sigma_P = 5$ dB

---

## 6. 안테나 모델링 (Section 7.3)

### BS 안테나 어레이 구조

$$
(M_g, N_g, M, N, P, d_H, d_V, d_{g,H}, d_{g,V})
$$

- $M_g \times N_g$ = 패널 수 (수직 × 수평)
- $M \times N$ = 패널당 엘리먼트 수 (수직 × 수평, 동일 편파)
- $P$ = 편파 수 (1 또는 2)
- $d_H$, $d_V$ = 엘리먼트 간격 [$\lambda$]
- $d_{g,H}$, $d_{g,V}$ = 패널 간격 [$\lambda$]

**엘리먼트 위치 계산** (Initiallization.cpp):

$$
d_{tx}[m][n][p][m_g][n_g].y = n \cdot d_H + n_g \cdot d_{g,H} \quad \text{(수평)}
$$

$$
d_{tx}[m][n][p][m_g][n_g].z = m \cdot d_V + m_g \cdot d_{g,V} \quad \text{(수직)}
$$

### 단일 엘리먼트 방사 패턴 (Table 7.3-1)

**수직 컷:** ($\theta_{3dB} = 65°$, GCS 기준 $\theta$는 zenith 각)

$$
A_{E,V}(\theta) = -\min\!\left\{ 12 \left(\frac{\theta - 90°}{\theta_{3dB}}\right)^2,\; SLA_V \right\}, \quad SLA_V = 30 \text{ dB}
$$

**수평 컷:** ($\varphi_{3dB} = 65°$)

$$
A_{E,H}(\varphi) = -\min\!\left\{ 12 \left(\frac{\varphi}{\varphi_{3dB}}\right)^2,\; A_{\max} \right\}, \quad A_{\max} = 30 \text{ dB}
$$

**3D 패턴:**

$$
A_E(\theta, \varphi) = -\min\!\left\{ -\left(A_{E,V}(\theta) + A_{E,H}(\varphi)\right),\; A_{\max} \right\}
$$

$$
G_{E,\max} = 8 \text{ dBi}
$$

### Antenna Port Mapping (Beamtilt)

$$
w_m = \frac{1}{\sqrt{M}} \exp\!\left(-j \frac{2\pi}{\lambda} (m-1) d_V \cos\theta_{etilt}\right) \tag{7.3-1}
$$

### Polarized Antenna Modelling

**Model-2** (calibration에서 사용):

$$
F_{\theta'}(\theta', \varphi') = \sqrt{A'(\theta', \varphi')} \cdot \cos(\zeta) \tag{7.3-4}
$$

$$
F_{\varphi'}(\theta', \varphi') = \sqrt{A'(\theta', \varphi')} \cdot \sin(\zeta) \tag{7.3-5}
$$

- $\zeta = 0°$: 수직 편파
- $\zeta = \pm 45°$: X-pol 쌍

**Model-1** (더 정확):

$$
\begin{bmatrix}
F_{\theta'} \\
F_{\varphi'}
\end{bmatrix}
=
\begin{bmatrix}
+\cos\psi & -\sin\psi \\
+\sin\psi & +\cos\psi
\end{bmatrix}
\begin{bmatrix}
F_{\theta''} \\
F_{\varphi''}
\end{bmatrix} \tag{7.3-3}
$$

### UT 안테나 (Handheld)

Directional 패턴 (Table 7.3-2):

$$
\theta_{3dB} = 90°, \quad SLA_V = 25 \text{ dB}, \quad \varphi_{3dB} = 90°, \quad A_{\max} = 25 \text{ dB}, \quad G_{E,\max} = 5 \text{ dBi}
$$

**Note:** Full calibration (Table 7.8-2)에서는 UT 안테나 패턴 = **Isotropic** ($G_{E,\max} = 0$ dBi)

---

## 7. LSP 파라미터 테이블 (Table 7.5-6)

### Part 1: UMi-Street Canyon & UMa

$f_c$: GHz, $d_{2D}$: km

#### Delay Spread (DS)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| $\mu_{\lg DS}$ | $-0.18\log_{10}(1+f_c)-7.28$ | $-0.22\log_{10}(1+f_c)-6.87$ | $-6.62$ | $-7.067-0.0794\log_{10}(f_c)$ | $-6.47-0.134\log_{10}(f_c)$ | $-6.62$ |
| $\sigma_{\lg DS}$ | $0.39$ | $0.19\log_{10}(1+f_c)+0.22$ | $0.32$ | $0.57+0.026\log_{10}(f_c)$ | $0.39$ | $0.32$ |

#### AOD Spread (ASD)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| $\mu_{\lg ASD}$ | $-0.05\log_{10}(1+f_c)+1.21$ | $-0.24\log_{10}(1+f_c)+1.54$ | $1.25$ | $0.92$ | $1.09$ | $0.58$ |
| $\sigma_{\lg ASD}$ | $0.08\log_{10}(1+f_c)+0.29$ | $0.10\log_{10}(1+f_c)+0.33$ | $0.42$ | $0.31$ | $0.44$ | $0.7$ |

#### AOA Spread (ASA)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| $\mu_{\lg ASA}$ | $-0.07\log_{10}(1+f_c)+1.66$ | $-0.07\log_{10}(1+f_c)+1.76$ | $1.76$ | $1.76$ | $2.04-0.25\log_{10}(f_c)$ | $1.76$ |
| $\sigma_{\lg ASA}$ | $0.021\log_{10}(1+f_c)+0.26$ | $0.05\log_{10}(1+f_c)+0.27$ | $0.16$ | $0.19$ | $0.17-0.03\log_{10}(f_c)$ | $0.16$ |

#### ZOA Spread (ZSA)

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| $\mu_{\lg ZSA}$ | $-0.11\log_{10}(1+f_c)+0.81$ | $-0.03\log_{10}(1+f_c)+0.92$ | $1.01$ | $0.96$ | $-0.2856\log_{10}(f_c)+1.445$ | $1.01$ |
| $\sigma_{\lg ZSA}$ | $-0.03\log_{10}(1+f_c)+0.29$ | $-0.05\log_{10}(1+f_c)+0.35$ | $0.43$ | $0.15$ | $0.17$ | $0.43$ |

#### K-factor, XPR, Clusters

| 파라미터 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| $\mu_K$ [dB] | 9 | N/A | N/A | 9 | N/A | N/A |
| $\sigma_K$ [dB] | 5 | N/A | N/A | 3.5 | N/A | N/A |
| $\mu_{\text{XPR}}$ [dB] | 9 | 8 | 9 | 8 | 7 | 9 |
| $\sigma_{\text{XPR}}$ [dB] | 3 | 3 | 5 | 4 | 3 | 5 |
| $N_{\text{clusters}}$ | 12 | 19 | 12 | 12 | 20 | 12 |
| $M_{\text{rays}}$ | 20 | 20 | 20 | 20 | 20 | 20 |
| $r_\tau$ | 3 | 2.1 | 2.2 | 2.5 | 2.3 | 2.2 |
| $c_{DS}$ [ns] | 5 | 11 | 11 | $\max(0.25,\; 6.5622-3.4084\log_{10}(f_c))$ | $\max(0.25,\; 6.5622-3.4084\log_{10}(f_c))$ | 11 |
| $c_{ASD}$ [°] | 3 | 10 | 5 | 3.58 | 1.8 | 1.8 |
| $c_{ASA}$ [°] | 17 | 22 | 8 | 11 | 15 | 8 |
| $c_{ZSA}$ [°] | 7 | 7 | 3 | 7 | 7 | 3 |
| $\zeta$ [dB] | 3 | 3 | 4 | 3 | 3 | 4 |

> **NOTE 8:** UMa $c_{ASD}$ (LOS=3.58, NLOS=1.8, O2I=1.8)는 Release 19에서 업데이트됨.

#### Cross-Correlations

| 쌍 | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| ASD-DS | 0.5 | 0 | 0.4 | 0.4 | 0.4 | 0.4 |
| ASA-DS | 0.8 | 0.4 | 0.4 | 0.8 | 0.6 | 0.4 |
| ASA-SF | -0.4 | -0.4 | 0 | -0.5 | 0 | 0 |
| ASD-SF | -0.5 | 0 | 0.2 | -0.5 | -0.6 | 0.2 |
| DS-SF | -0.4 | -0.7 | -0.5 | -0.4 | -0.4 | -0.5 |
| ASD-ASA | 0.4 | 0 | 0 | 0 | 0.4 | 0 |
| ASD-K | -0.2 | N/A | N/A | 0 | N/A | N/A |
| ASA-K | -0.3 | N/A | N/A | -0.2 | N/A | N/A |
| DS-K | -0.7 | N/A | N/A | -0.4 | N/A | N/A |
| SF-K | 0.5 | N/A | N/A | 0 | N/A | N/A |
| ZSD-SF | 0 | 0 | 0 | 0 | 0 | 0 |
| ZSA-SF | 0 | 0 | 0 | -0.8 | -0.4 | 0 |
| ZSD-K | 0 | N/A | N/A | 0 | N/A | N/A |
| ZSA-K | 0 | N/A | N/A | 0 | N/A | N/A |
| ZSD-DS | 0 | -0.5 | -0.6 | -0.2 | -0.5 | -0.6 |
| ZSA-DS | 0.2 | 0 | -0.2 | 0 | 0 | -0.2 |
| ZSD-ASD | 0.5 | 0.5 | -0.2 | 0.5 | 0.5 | -0.2 |
| ZSA-ASD | 0.3 | 0.5 | 0 | 0 | -0.1 | 0 |
| ZSD-ASA | 0 | 0 | 0 | -0.3 | 0 | 0 |
| ZSA-ASA | 0 | 0.2 | 0.5 | 0.4 | 0 | 0.5 |
| ZSD-ZSA | 0 | 0 | 0.5 | 0 | 0 | 0.5 |

#### Correlation Distance [m]

| LSP | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I |
|---|---|---|---|---|---|---|
| DS | 7 | 10 | 10 | 30 | 40 | 10 |
| ASD | 8 | 10 | 11 | 18 | 50 | 11 |
| ASA | 8 | 9 | 17 | 15 | 50 | 17 |
| SF | 10 | 13 | 7 | 37 | 50 | 7 |
| K | 15 | N/A | N/A | 12 | N/A | N/A |
| ZSA | 12 | 10 | 25 | 15 | 50 | 25 |
| ZSD | 12 | 10 | 25 | 15 | 50 | 25 |

### Part 2: RMa (up to 7 GHz) & Indoor-Office

#### Delay/Angle Spread

| 파라미터 | RMa LOS | RMa NLOS | RMa O2I | InH LOS | InH NLOS |
|---|---|---|---|---|---|
| $\mu_{\lg DS}$ | $-7.49$ | $-7.43$ | $-7.47$ | $-0.01\log_{10}(1+f_c)-7.692$ | $-0.28\log_{10}(1+f_c)-7.173$ |
| $\sigma_{\lg DS}$ | $0.55$ | $0.48$ | $0.24$ | $0.18$ | $0.10\log_{10}(1+f_c)+0.055$ |
| $\mu_{\lg ASD}$ | $0.90$ | $0.95$ | $0.67$ | $1.60$ | $1.62$ |
| $\sigma_{\lg ASD}$ | $0.38$ | $0.45$ | $0.18$ | $0.18$ | $0.25$ |
| $\mu_{\lg ASA}$ | $1.52$ | $1.52$ | $1.66$ | $-0.19\log_{10}(1+f_c)+1.781$ | $-0.11\log_{10}(1+f_c)+1.863$ |
| $\sigma_{\lg ASA}$ | $0.24$ | $0.13$ | $0.21$ | $0.12\log_{10}(1+f_c)+0.119$ | $0.12\log_{10}(1+f_c)+0.059$ |
| $\mu_{\lg ZSA}$ | $0.47$ | $0.58$ | $0.93$ | $-0.26\log_{10}(1+f_c)+1.44$ | $-0.15\log_{10}(1+f_c)+1.387$ |
| $\sigma_{\lg ZSA}$ | $0.40$ | $0.37$ | $0.22$ | $-0.04\log_{10}(1+f_c)+0.264$ | $-0.09\log_{10}(1+f_c)+0.746$ |

#### K-factor, SF, XPR, Clusters

| 파라미터 | RMa LOS | RMa NLOS | RMa O2I | InH LOS | InH NLOS |
|---|---|---|---|---|---|
| $\sigma_{SF}$ [dB] | Table 7.4.1-1 | Table 7.4.1-1 | 8 | Table 7.4.1-1 | Table 7.4.1-1 |
| $\mu_K$ [dB] | 7 | N/A | N/A | 7 | N/A |
| $\sigma_K$ [dB] | 4 | N/A | N/A | 4 | N/A |
| $\mu_{\text{XPR}}$ [dB] | 12 | 7 | 7 | 11 | 10 |
| $\sigma_{\text{XPR}}$ [dB] | 4 | 3 | 3 | 4 | 4 |
| $N_{\text{clusters}}$ | 11 | 10 | 10 | 15 | 19 |
| $M_{\text{rays}}$ | 20 | 20 | 20 | 20 | 20 |
| $r_\tau$ | 3.8 | 1.7 | 1.7 | 3.6 | 3 |
| $c_{DS}$ [ns] | N/A | N/A | N/A | N/A | N/A |
| $c_{ASD}$ [°] | 2 | 2 | 2 | 5 | 5 |
| $c_{ASA}$ [°] | 3 | 3 | 3 | 8 | 11 |
| $c_{ZSA}$ [°] | 3 | 3 | 3 | 9 | 9 |
| $\zeta$ [dB] | 3 | 3 | 3 | 6 | 3 |

#### Cross-Correlations (RMa / InH)

| 쌍 | RMa LOS | RMa NLOS | RMa O2I | InH LOS | InH NLOS |
|---|---|---|---|---|---|
| ASD-DS | 0 | -0.4 | 0 | 0.6 | 0.4 |
| ASA-DS | 0 | 0 | 0 | 0.8 | 0 |
| ASA-SF | 0 | 0 | 0 | -0.5 | -0.4 |
| ASD-SF | 0 | 0.6 | 0 | -0.4 | 0 |
| DS-SF | -0.5 | -0.5 | 0 | -0.8 | -0.5 |
| ASD-ASA | 0 | 0 | -0.7 | 0.4 | 0 |
| ASD-K | 0 | N/A | N/A | 0 | N/A |
| ASA-K | 0 | N/A | N/A | 0 | N/A |
| DS-K | 0 | N/A | N/A | -0.5 | N/A |
| SF-K | 0 | N/A | N/A | 0.5 | N/A |
| ZSD-SF | 0.01 | -0.04 | 0 | 0.2 | 0 |
| ZSA-SF | -0.17 | -0.25 | 0 | 0.3 | 0 |
| ZSD-K | 0 | N/A | N/A | 0 | N/A |
| ZSA-K | -0.02 | N/A | N/A | 0.1 | N/A |
| ZSD-DS | -0.05 | -0.10 | 0 | 0.1 | -0.27 |
| ZSA-DS | 0.27 | -0.40 | 0 | 0.2 | -0.06 |
| ZSD-ASD | 0.73 | 0.42 | 0.66 | 0.5 | 0.35 |
| ZSA-ASD | -0.14 | -0.27 | 0.47 | 0 | 0.23 |
| ZSD-ASA | -0.20 | -0.18 | -0.55 | 0 | -0.08 |
| ZSA-ASA | 0.24 | 0.26 | -0.22 | 0.5 | 0.43 |
| ZSD-ZSA | -0.07 | -0.27 | 0 | 0 | 0.42 |

#### Correlation Distance [m] (RMa / InH)

| LSP | RMa LOS | RMa NLOS | RMa O2I | InH LOS | InH NLOS |
|---|---|---|---|---|---|
| DS | 50 | 36 | 36 | 8 | 5 |
| ASD | 25 | 30 | 30 | 7 | 3 |
| ASA | 35 | 40 | 40 | 5 | 3 |
| SF | 37 | 120 | 120 | 10 | 6 |
| K | 40 | N/A | N/A | 4 | N/A |
| ZSA | 15 | 50 | 50 | 4 | 4 |
| ZSD | 15 | 50 | 50 | 4 | 4 |

> **NOTE:** InH에서 $f_c < 6$ GHz이면 $f_c = 6$으로 치환하여 LSP 결정

### Part 3: InF (Indoor Factory)

> $V$/$S$: 홀 체적(m³)/전체 표면적(m², 벽+바닥+천장)

| 파라미터 | InF LOS | InF NLOS |
|---|---|---|
| $\mu_{\lg DS}$ | $\log_{10}(26(V/S)+14)-9.35$ | $\log_{10}(30(V/S)+32)-9.44$ |
| $\sigma_{\lg DS}$ | $0.15$ | $0.19$ |
| $\mu_{\lg ASD}$ | $1.56$ | $1.57$ |
| $\sigma_{\lg ASD}$ | $0.25$ | $0.2$ |
| $\mu_{\lg ASA}$ | $-0.18\log_{10}(1+f_c)+1.78$ | $1.72$ |
| $\sigma_{\lg ASA}$ | $0.12\log_{10}(1+f_c)+0.2$ | $0.3$ |
| $\mu_{\lg ZSA}$ | $-0.2\log_{10}(1+f_c)+1.5$ | $-0.13\log_{10}(1+f_c)+1.45$ |
| $\sigma_{\lg ZSA}$ | $0.35$ | $0.45$ |
| $\sigma_{SF}$ [dB] | Table 7.4.1-1 | Table 7.4.1-1 |
| $\mu_K$ [dB] | 7 | N/A |
| $\sigma_K$ [dB] | 8 | N/A |
| $\mu_{\text{XPR}}$ [dB] | 12 | 11 |
| $\sigma_{\text{XPR}}$ [dB] | 6 | 6 |
| $N_{\text{clusters}}$ | 25 | 25 |
| $M_{\text{rays}}$ | 20 | 20 |
| $r_\tau$ | 2.7 | 3 |
| $c_{DS}$ [ns] | N/A | N/A |
| $c_{ASD}$ [°] | 5 | 5 |
| $c_{ASA}$ [°] | 8 | 8 |
| $c_{ZSA}$ [°] | 9 | 9 |
| $\zeta$ [dB] | 4 | 3 |

**InF Cross-Correlations:** 모든 쌍이 **0** (상관 없음). 단, ASD-K=-0.5, DS-K=-0.7 (LOS).

**InF Correlation Distance:** 모두 **10 m**.

### Part 4: SMa (Suburban Macro, up to 24 GHz, V19 업데이트)

> **New in Rel-19.** SMa 시나리오의 LSP 파라미터.

| 파라미터 | SMa LOS | SMa NLOS | SMa O2I |
|---|---|---|---|
| $\mu_{\lg DS}$ | $-7.42$ | $-7.20$ | $-7.20$ |
| $\sigma_{\lg DS}$ | $0.60$ | $0.58$ | $0.58$ |
| $\mu_{\lg ASD}$ | $0.48$ | $0.51$ | $0.51$ |
| $\sigma_{\lg ASD}$ | $0.27$ | $0.33$ | $0.33$ |
| $\mu_{\lg ASA}$ | $1.43$ | $1.63$ | $1.63$ |
| $\sigma_{\lg ASA}$ | $0.12$ | $0.26$ | $0.26$ |
| $\mu_{\lg ZSA}$ | $1.18$ | $1.16$ | $1.16$ |
| $\sigma_{\lg ZSA}$ | $0.05$ | $0.14$ | $0.14$ |
| $\sigma_{SF}$ [dB] | Table 7.4.1-1 | Table 7.4.1-1 | 8 |
| $\mu_K$ [dB] | 9 | N/A | N/A |
| $\sigma_K$ [dB] | 7 | N/A | N/A |
| $\mu_{\text{XPR}}$ [dB] | 8 | 4 | 4 |
| $\sigma_{\text{XPR}}$ [dB] | 4 | 3 | 3 |
| $N_{\text{clusters}}$ | 15 | 14 | 14 |
| $M_{\text{rays}}$ | 20 | 20 | 20 |
| $r_\tau$ | 2.4 | 1.5 | 1.5 |
| $c_{DS}$ [ns] | $\max(0.25,\; 6.5622-3.4084\log_{10}(f_c))$ | $\max(0.25,\; 6.5622-3.4084\log_{10}(f_c))$ | $\max(0.25,\; 6.5622-3.4084\log_{10}(f_c))$ |
| $c_{ASD}$ [°] | 2.08 | 1.33 | 1.33 |
| $c_{ASA}$ [°] | 5 | 10 | 10 |
| $c_{ZSA}$ [°] | 7 | 7 | 7 |
| $\zeta$ [dB] | 3 | 3 | 3 |

#### SMa Cross-Correlations

| 쌍 | SMa LOS | SMa NLOS | SMa O2I |
|---|---|---|---|
| ASD-DS | 0.4 | 0 | 0 |
| ASA-DS | 0.8 | 0.7 | 0.7 |
| ASA-SF | -0.5 | 0 | 0 |
| ASD-SF | -0.5 | -0.4 | -0.4 |
| DS-SF | -0.4 | -0.13 | -0.13 |
| ASD-ASA | 0 | 0 | 0 |
| ASD-K | 0 | N/A | N/A |
| ASA-K | -0.2 | N/A | N/A |
| DS-K | -0.4 | N/A | N/A |
| SF-K | 0 | N/A | N/A |
| ZSD-SF | 0 | 0 | 0 |
| ZSA-SF | -0.8 | -0.47 | -0.47 |
| ZSD-K | 0 | N/A | N/A |
| ZSA-K | 0 | N/A | N/A |
| ZSD-DS | -0.2 | -0.5 | -0.5 |
| ZSA-DS | 0 | 0 | 0 |
| ZSD-ASD | 0.5 | 0.5 | 0.5 |
| ZSA-ASD | 0 | -0.1 | -0.1 |
| ZSD-ASA | -0.3 | 0 | 0 |
| ZSA-ASA | 0.4 | 0.36 | 0.36 |
| ZSD-ZSA | 0 | 0 | 0 |

#### SMa Correlation Distance [m]

| LSP | SMa LOS | SMa NLOS | SMa O2I |
|---|---|---|---|
| DS | 6 | 40 | 40 |
| ASD | 15 | 30 | 30 |
| ASA | 20 | 30 | 30 |
| SF | 40 | 50 | 50 |
| K | 10 | N/A | N/A |
| ZSA | 15 | 50 | 50 |
| ZSD | 15 | 50 | 50 |

> **주파수 치환 규칙:**
> - UMa: $f_c < 6$ GHz → $f_c = 6$으로 치환
> - UMi: $f_c < 2$ GHz → $f_c = 2$로 치환
> - InH: $f_c < 6$ GHz → $f_c = 6$으로 치환

---

## 8. ZSD/ZOD 오프셋 (Table 7.5-7~12)

### UMa (Table 7.5-7)

| | LOS | NLOS |
|---|---|---|
| $\mu_{\lg ZSD}$ | $\max\!\left[-0.5,\; -2.1(d_{2D}/1000) - 0.01(h_{UT}-1.5) + 0.75\right]$ | $\max\!\left[-0.5,\; -2.1(d_{2D}/1000) - 0.01(h_{UT}-1.5) + 0.9\right]$ |
| $\sigma_{\lg ZSD}$ | 0.40 | 0.49 |
| $\mu_{\text{offset,ZOD}}$ | 0 | $e(f_c) - 10^{a(f_c)\log_{10}(\max(b(f_c), d_{2D})) + c(f_c) - 0.07(h_{UT}-1.5)}$ |

NLOS ZOD offset 계수:

$$
a(f_c) = 0.208\log_{10}(f_c) - 0.782
$$

$$
b(f_c) = 25
$$

$$
c(f_c) = -0.13\log_{10}(f_c) + 2.03
$$

$$
e(f_c) = 7.66\log_{10}(f_c) - 5.96
$$

### UMi (Table 7.5-8)

| | LOS | NLOS |
|---|---|---|
| $\mu_{\lg ZSD}$ | $\max\!\left[-0.21,\; -14.8(d_{2D}/1000) + 0.01|h_{UT}-h_{BS}| + 0.83\right]$ | $\max\!\left[-0.5,\; -3.1(d_{2D}/1000) + 0.01\max(h_{UT}-h_{BS},0) + 0.2\right]$ |
| $\sigma_{\lg ZSD}$ | 0.35 | 0.35 |
| $\mu_{\text{offset,ZOD}}$ | 0 | $-10^{-1.5\log_{10}(\max(10, d_{2D})) + 3.3}$ |

### RMa (Table 7.5-9)

| | LOS | NLOS / O2I |
|---|---|---|
| $\mu_{\lg ZSD}$ | $\max\!\left[-1,\; -0.17(d_{2D}/1000) - 0.01(h_{UT}-1.5) + 0.22\right]$ | $\max\!\left[-1,\; -0.19(d_{2D}/1000) - 0.01(h_{UT}-1.5) + 0.28\right]$ |
| $\sigma_{\lg ZSD}$ | 0.34 | 0.30 |
| $\mu_{\text{offset,ZOD}}$ | 0 | $\arctan\!\left(\frac{35-3.5}{d_{2D}}\right) - \arctan\!\left(\frac{35-1.5}{d_{2D}}\right)$ |

### InH (Table 7.5-10)

| | LOS | NLOS |
|---|---|---|
| $\mu_{\lg ZSD}$ | $-1.43\log_{10}(1+f_c) + 2.228$ | $1.08$ |
| $\sigma_{\lg ZSD}$ | $0.13\log_{10}(1+f_c) + 0.30$ | $0.36$ |
| $\mu_{\text{offset,ZOD}}$ | 0 | 0 |

### InF (Table 7.5-11)

| | LOS | NLOS |
|---|---|---|
| $\mu_{\lg ZSD}$ | $1.35$ | $1.2$ |
| $\sigma_{\lg ZSD}$ | $0.35$ | $0.55$ |
| $\mu_{\text{offset,ZOD}}$ | 0 | 0 |

### SMa (Table 7.5-12, V19 업데이트)

| | LOS | NLOS | O2I |
|---|---|---|---|
| $\mu_{\lg ZSD}$ | 0.14 | 0.14 | 0.14 |
| $\sigma_{\lg ZSD}$ | 0.16 | 0.16 | 0.16 |
| $\mu_{\text{offset,ZOD}}$ | 0 | 3.5 | 3.5 |

> **ZSD/ZOD Tables 공통 Notes:**
> - $f_c$: 중심 주파수 [GHz], $d_{2D}$: BS-UT 거리 [m], $h_{BS}$/$h_{UT}$: 안테나 높이 [m]
> - 표기법: $\mu_{\lg X} = \text{mean}\{\log_{10}(X)\}$, $\sigma_{\lg X} = \text{std}\{\log_{10}(X)\}$
> - **$f_c < 6$ GHz일 때 Table 7.5-7과 7.5-10의 frequency-dependent 파라미터는 $f_c = 6$으로 치환**
> - **O2I 링크의 ZSD 파라미터:** outdoor 링크 부분의 LOS 조건에 따라 Table 7.5-7/8의 outdoor 파라미터를 사용

---

## 9. Scaling Factor 테이블

(Section 7 Step 7에 통합됨. Table 7.5-2, 7.5-3, 7.5-4, 7.5-5 참조.)

---

## 10. Calibration 파라미터 (Section 7.8)

### Table 7.8-1: Large Scale Calibration

| Parameter | Values |
|---|---|
| Scenarios | UMa, UMi-Street Canyon, Indoor-office (open office) |
| Sectorization | 3 sectors: 30°, 150°, 270° |
| BS antenna | $M_g=N_g=1$, $(M,N,P)=(10,1,1)$, $d_V=0.5\lambda$ |
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
| **BS Config 1** | $M=4, N=4, P=2, M_g=1, N_g=2$, $d_H=d_V=0.5\lambda$, $d_{H,g}=d_{V,g}=2.5\lambda$ → metrics 1,2,3 |
| BS Config 1 port | 16 elem/pol/panel → 1 CRS port, panning (0,0)° |
| **BS Config 2** | $M_g=N_g=1$, $M=N=2$, $P=1$ → metrics 1,2,4 |
| BS Config 2 port | Each element → 1 CRS port |
| UT antenna | $M_g=N_g=1$, $M=N=1$, $P=2$ |
| UT attachment | RSRP (formula) from CRS port 0 |
| Polarized antenna | Model-2 (TR36.873) |
| UT orientation | $\Omega_\alpha \sim U[0,360°]$, $\Omega_\beta = 90°$, $\Omega_\gamma = 0°$ |
| UT antenna pattern | Isotropic |
| Polarization | $P=2$ → BS: X-pol ($\pm 45°$), UT: X-pol ($0°/+90°$) |
| Metrics | 1) Coupling loss, 2) Wideband SIR, 3) DS/AS CDFs, 4) PRB SVD CDFs |

---

## 11. 추가 모델링 컴포넌트 (Section 7.6)

> 이 섹션은 고급 시뮬레이션(대규모 어레이, 넓은 대역폭, 53-67 GHz 산소 흡수, 공간 일관성, 이동성, blockage 등)을 위한 추가 모델링이다.

### 7.6.1 산소 흡수 (Oxygen Absorption)

Step 11에서 생성된 클러스터 응답에 추가 손실을 적용한다:

$$
OL_n(f_c) = \frac{\alpha(f_c)}{1000} \cdot (d_{3D} + c \cdot \tau_n - c \cdot \tau_\Delta) \;\text{[dB]} \tag{7.6-1}
$$

- $\alpha(f_c)$: 주파수 의존 산소 손실 [dB/km] (Table 7.6.1-1)
- $d_{3D}$: BS-UT 거리 [m], $c$: 광속 [m/s]
- $\tau_n$: $n$번째 클러스터 delay [s]
- $\tau_\Delta = 0$ (LOS), $\min(\tau_n)$ (NLOS)

**Table 7.6.1-1: 주파수별 산소 흡수 손실 $\alpha(f)$ [dB/km]**

| $f$ [GHz] | 0-52 | 53 | 54 | 55 | 56 | 57 | 58 | 59 | 60 | 61 | 62 | 63 | 64 | 65 | 66-100 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $\alpha$ [dB/km] | 0 | 1 | 2.2 | 4 | 6.6 | 9.7 | 12.6 | 14.6 | 15 | 14.6 | 14.1 | 10.7 | 6.7 | 3.5 | 1 |

> 테이블에 없는 주파수는 양쪽 인접 값에서 **선형 보간**한다.

### 7.6.2 대역폭/어레이 확장 (Large Bandwidth and Large Array)

대역폭 $B > c/D$ Hz ($D$: 최대 안테나 aperture)인 경우, 각 ray에 대해 개별 TOA를 모델링한다.

### 7.6.3 공간 일관성 (Spatial Consistency)

LSP 및 SSP에 공간 일관성을 적용하여 인접 UT 간 채널 특성의 연속적 변화를 모델링한다.
- LSP 상관 거리: Table 7.5-6의 correlation distance 사용
- 클러스터 delay/power/angle의 공간 상관: Table 7.6.3.1-2의 상관 거리 적용

**Table 7.6.3.1-2: Cluster 파라미터 공간 상관 거리 [m]**

| 파라미터 | RMa LOS | RMa NLOS | RMa O2I | UMi LOS | UMi NLOS | UMi O2I | UMa LOS | UMa NLOS | UMa O2I | InH LOS | InH NLOS | InF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Delay | 50 | 36 | 36 | 12 | 15 | 15 | 7 | 11 | 11 | 10 | 10 | 10 |
| AOD/AOA | 50 | 36 | 36 | 12 | 15 | 15 | 7 | 11 | 11 | 10 | 10 | 10 |
| ZOD/ZOA | 50 | 36 | 36 | 12 | 15 | 15 | 7 | 11 | 11 | 10 | 10 | 10 |
| Shadow | 50 | 36 | 36 | 12 | 15 | 15 | 7 | 11 | 11 | 10 | 10 | 10 |

### 7.6.4 Blockage

Step 11의 채널 계수에 blockage 모델을 적용한다. Blocker는 사각형 스크린(화면)으로 모델링하며, 각 클러스터/ray에 대해 attenuation을 계산한다.

### 7.6.5 Shadow Fading 자기상관

거리 $x$에 대한 SF 자기상관:

$$
R(x) = e^{-x/d_{\text{cor}}} \tag{7.4-5}
$$

여기서 $d_{\text{cor}}$은 상관 거리 (Table 7.5-6의 SF correlation distance).

---

## 12. RSRP 계산 공식 (TR 36.873)

> **Source:** 3GPP TR 36.873 V12.7.0 (2017-12), Release 12, Section 8.1

Channel model calibration 및 baseline performance evaluation을 위한 RSRP 계산 공식이다. TX 안테나 port $p$에 대해 UE attachment에 필요한 RSRP를 계산한다 (NLOS UE에는 $K_R = 0$ 적용). TX power는 CRS transmitted power per RE이다. $w_s$ ($s=1,\ldots,S$)는 port $p$의 virtualization을 위한 complex weight vector이며, $U$는 수신 안테나 엘리먼트 수이다.

$$
RSRP_p = PL \cdot SF \cdot \sum_{u=1}^{U} \left( |\alpha_{0,u,p}|^2 + \sum_{n=1}^{N} \sum_{m=1}^{M} |\alpha_{n,m,u,p}|^2 \right) \cdot \frac{TX_{\text{power}}}{U} \tag{8.1-1}
$$

여기서 NLOS path ($n=1,\ldots,N$)에 대해:

$$
\alpha_{n,m,u,p} = \sqrt{\frac{P_n}{M(K_R+1)}}
\begin{bmatrix}
F_{rx,u,\theta}(\theta_{n,m,\text{ZOA}}, \phi_{n,m,\text{AOA}}) \\
F_{rx,u,\varphi}(\theta_{n,m,\text{ZOA}}, \phi_{n,m,\text{AOA}})
\end{bmatrix}^T
\begin{bmatrix}
e^{j\Phi^{\theta\theta}_{n,m}} & \sqrt{\kappa_{n,m}^{-1}} e^{j\Phi^{\theta\varphi}_{n,m}} \\
\sqrt{\kappa_{n,m}^{-1}} e^{j\Phi^{\varphi\theta}_{n,m}} & e^{j\Phi^{\varphi\varphi}_{n,m}}
\end{bmatrix}
\begin{bmatrix}
F_{tx,p,\theta}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}}) \\
F_{tx,p,\varphi}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}})
\end{bmatrix} \tag{8.1-2}
$$

LOS path에 대해:

$$
\alpha_{0,u,p} = \sqrt{\frac{K_R}{K_R+1}}
\begin{bmatrix}
F_{rx,u,\theta}(\theta_{\text{LOS,ZOA}}, \phi_{\text{LOS,AOA}}) \\
F_{rx,u,\varphi}(\theta_{\text{LOS,ZOA}}, \phi_{\text{LOS,AOA}})
\end{bmatrix}^T
\begin{bmatrix}
e^{j\Phi_{\text{LOS}}} & 0 \\
0 & -e^{j\Phi_{\text{LOS}}}
\end{bmatrix}
\begin{bmatrix}
F_{tx,p,\theta}(\theta_{\text{LOS,ZOD}}, \phi_{\text{LOS,AOD}}) \\
F_{tx,p,\varphi}(\theta_{\text{LOS,ZOD}}, \phi_{\text{LOS,AOD}})
\end{bmatrix} \tag{8.1-3}
$$

여기서 port-level field pattern:

$$
F_{tx,p,\theta}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}}) = \sum_{s=1}^{S} w_s \, e^{j \frac{2\pi}{\lambda_0} \hat{r}_{tx,n,m}^T \cdot \bar{d}_{tx,s}} \, F_{tx,s,\theta}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}}) \tag{8.1-4}
$$

$$
F_{tx,p,\varphi}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}}) = \sum_{s=1}^{S} w_s \, e^{j \frac{2\pi}{\lambda_0} \hat{r}_{tx,n,m}^T \cdot \bar{d}_{tx,s}} \, F_{tx,s,\varphi}(\theta_{n,m,\text{ZOD}}, \phi_{n,m,\text{AOD}}) \tag{8.1-5}
$$

### 변수 설명

| 기호 | 설명 |
|------|------|
| $RSRP_p$ | TX 안테나 port $p$에 대한 RSRP |
| $PL$ | Pathloss (선형값) |
| $SF$ | Shadow Fading (선형값) |
| $U$ | RX 안테나 엘리먼트 수 |
| $N$ | 클러스터 수 |
| $M$ | 클러스터당 레이 수 |
| $P_n$ | 클러스터 $n$의 파워 (K-factor 없이 정규화, $\sum P_n = 1$) |
| $K_R$ | Ricean K-factor (선형값, NLOS일 경우 0) |
| $w_s$ | Port $p$의 virtualization weight (엘리먼트 $s$에 대한 복소 가중치) |
| $S$ | Port $p$에 매핑된 TX 엘리먼트 수 |
| $F_{rx,u,\theta/\varphi}$ | RX 엘리먼트 $u$의 field pattern ($\theta$, $\varphi$ 성분) |
| $F_{tx,s,\theta/\varphi}$ | TX 엘리먼트 $s$의 element-level field pattern |
| $F_{tx,p,\theta/\varphi}$ | TX port $p$의 port-level field pattern (식 8.1-4, 8.1-5) |
| $\kappa_{n,m}$ | Cross-polarization ratio (XPR) |
| $\Phi^{xx}_{n,m}$ | 랜덤 초기 위상 ($\theta\theta$, $\theta\varphi$, $\varphi\theta$, $\varphi\varphi$) |
| $\hat{r}_{tx,n,m}$ | 레이 $(n,m)$의 TX 출발 방향 단위벡터 |
| $\bar{d}_{tx,s}$ | TX 엘리먼트 $s$의 위치 벡터 |
| $\lambda_0$ | 반송파 파장 |

---

## 13. Annex A: Circular Angle Spread 계산

> **Source:** 3GPP TR 25.996 V19.0.0, Annex A

$N$개의 multi-path, 각각 $M$개의 sub-path를 가진 신호에 대한 composite signal의 conventional angle spread 계산:

$$
\sigma_{AS} = \sqrt{\frac{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} (\theta_{n,m,\mu})^2 \cdot P_{n,m}}{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} P_{n,m}}} \tag{A-1}
$$

여기서 $P_{n,m}$은 $n$번째 path의 $m$번째 subpath의 power이고, $\theta_{n,m,\mu}$는 다음과 같이 정의된다:

$$
\theta_{n,m,\mu} = \text{mod}\!\left(\theta_{n,m} - \mu_\theta + \pi,\; 2\pi\right) - \pi \tag{A-2}
$$

$\mu_\theta$는 다음과 같이 정의된다:

$$
\mu_\theta = \frac{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} \theta_{n,m} \cdot P_{n,m}}{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} P_{n,m}} \tag{A-3}
$$

$\theta_{n,m}$은 $n$번째 path의 $m$번째 subpath의 AoA (또는 AoD)이다.

Angle spread는 AoA의 선형 이동에 독립적이어야 한다. 즉, $\theta_{n,m}$을 $\theta_{n,m}+\Delta$로 대체해도 $\sigma_{AS}(\Delta)$는 일정해야 한다. 그러나 modulo $2\pi$ 연산의 모호성으로 인해 그렇지 않을 수 있으므로, angle spread는 모든 $\Delta$에 대한 $\sigma_{AS}(\Delta)$의 최솟값이어야 한다:

$$
\sigma_{AS} = \min_\Delta \sigma_{AS}(\Delta) = \sqrt{\frac{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} (\theta_{n,m,\mu}(\Delta))^2 \cdot P_{n,m}}{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} P_{n,m}}} \tag{A-4}
$$

여기서:

$$
\theta_{n,m,\mu}(\Delta) = \text{mod}\!\left(\theta_{n,m}(\Delta) - \mu_\theta(\Delta) + \pi,\; 2\pi\right) - \pi \tag{A-5}
$$

$$
\mu_\theta(\Delta) = \frac{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} \theta_{n,m}(\Delta) \cdot P_{n,m}}{\displaystyle\sum_{n=1}^{N}\sum_{m=1}^{M} P_{n,m}} \tag{A-6}
$$

그리고 $\theta_{n,m}(\Delta) = \text{mod}(\theta_{n,m} + \Delta + \pi,\; 2\pi) - \pi$.
