# Spatial Non-Stationarity (SNS) 구현 가이드

## 참조 문서

- **3GPP TR 38.901 Release 20** (`doc/38901-j20.docx`)
  - Section 7.6.14.1.3: BS-side Stochastic Based Model
  - Section 7.6.14.2: UT-side SNS
  - Section 7.8.3: SNS Calibration

## 관련 코드 파일

| 파일 | 역할 |
|------|------|
| `h/channel.h` | `ClusterVR` 구조체, `compute_sns_attenuation()` 인라인 함수 |
| `src/channel.cpp` | `Generate_VisibilityRegion()`, `GetChannelImpulseResponse()` |
| `src/Link.cpp` | `find_best_tx_beam()`, `Get_BS_antenna_field_pattern()`, `Get_RSRP()` |
| `h/Initiallization.h` | SNS 글로벌 파라미터 (`g_sns_*`) |

---

## 1. 데이터 구조

### ClusterVR (h/channel.h)

```cpp
struct ClusterVR {
    bool  limited;    // true = visibility-limited cluster
    Real  V_n;        // VR size factor [0,1] (Eq. 7.6-57)
    Real  a;          // VR width (wavelength 단위)
    Real  b;          // VR height (wavelength 단위)
    Real  x0, y0;     // reference corner (VR가 anchored된 어레이 코너)
    Real  xa, yb;     // VR 내부 경계 (anchor 반대쪽)
    Real  xA, yB;     // 어레이 먼 코너 (anchor 대각선 반대)
    Real  D_n;        // sqrt((xA-xa)² + (yB-yb)²) — roll-off 정규화 거리
};
```

### CHANNEL 클래스 멤버

```cpp
ClusterVR sns_vr[MAX_NUM_CLUSTERS];  // NLOS 클러스터별 VR (서브클러스터 포함)
ClusterVR sns_vr_los;                // LOS path VR (= sns_vr[0] 복사)
bool      sns_any_limited;           // 하나라도 limited 클러스터가 있으면 true
Real      sns_Pr_sns;                // per-UT visibility probability
```

---

## 2. VR 생성: Generate_VisibilityRegion()

호출 위치: `Set_SmallScaleParameter()` → `Set_POWER()` → **`Generate_VisibilityRegion()`** → `Set_SUBCLUSTER()`

### Step 1: Per-UT 가시성 확률

```
Pr_sns = clamp(N(μ_vis, σ²_vis), 0, 1)    (Eq. 7.6-54)
```

### Step 2: 어레이 크기 (wavelength 단위)

```
W = (N-1)·dH + (Ng-1)·(N·dH + dgH)    // 수평
H = (M-1)·dV + (Mg-1)·(M·dV + dgV)    // 수직
```

### Step 3: 클러스터 파워와 P_max

**중요: `powerForAngles[n]`을 사용 (K-factor 반영됨)**

LOS 조건에서 powerForAngles는 Eq. 7.5-8에 따라 K-factor가 적용됨:
- 첫 번째 클러스터: `powerForAngles[0] = power[0]/(K+1) + K/(K+1)`
- 나머지 클러스터: `powerForAngles[n] = power[n]/(K+1)`
- NLOS 조건: `powerForAngles[n] = power[n]` (K-factor 없음)

```
P_max = max(powerForAngles[n])    // 모든 클러스터 중 최대값
```

LOS 조건에서는 powerForAngles[0]에 K/(K+1)이 이미 포함되어 있으므로,
LOS path를 별도의 "additional cluster"로 분리하여 P_max에 추가할 필요 없음.

### Step 4: 각 클러스터의 visibility 판정

```
x_n ~ U(0,1)
if x_n < Pr_sns → cluster n은 visibility-limited
```

### Step 5: VR 크기 — V_n (Eq. 7.6-57)

```
V_n = A · exp(-(P_max_dB - P_n_dB) / R) + B + δ    (Eq. 7.6-57)
V_n = clamp(V_n, 0, 1)
```

여기서:
- `P_n_dB = 10·log10(powerForAngles[n])`
- `P_max_dB = 10·log10(P_max)`
- `δ ~ N(0, σ_δ)` — δ는 **exp 바깥**에 더해짐 (additive)
- A, B, R: Table 7.6.14.1.2-2 파라미터

**주의**: δ를 exp() 안에 넣으면 안 됨. `A·exp(...+δ)`가 아닌 `A·exp(...) + B + δ`가 올바름.

### Step 6: VR 크기와 위치

```
a ~ U(V_n·W, W)              // VR width
b = V_n·H·W / a              // VR height (면적 = V_n·W·H 보존)
if b > H: b = H, a 재조정
```

VR은 어레이의 4개 코너 중 하나에 랜덤 배치 (각 50% 확률로 좌/우, 상/하).

### Step 7: LOS path VR

```cpp
// LOS path는 cluster 0과 동일한 공간 방향 → 같은 VR 공유
if (LOS) sns_vr_los = sns_vr[0];
```

powerForAngles[0]에 이미 K/(K+1)이 포함되어 있으므로 cluster 0이 가장 큰 V_n을 받음.
LOS와 cluster 0은 동일한 방향이므로 별도의 VR 계산이 불필요함.

### Step 8: 서브클러스터 VR 상속

`Set_SUBCLUSTER()`에서 2개의 strongest cluster가 3개의 서브클러스터로 분할됨.
서브클러스터(인덱스 `num_path ~ num_path+3`)는 부모 클러스터의 VR을 상속:

```cpp
// Set_SUBCLUSTER() 내:
sns_vr[num_path]     = sns_vr[strongest_power_idx];   // subcluster 2 of 1st strongest
sns_vr[num_path + 1] = sns_vr[strongest_power_idx];   // subcluster 3 of 1st strongest
sns_vr[num_path + 2] = sns_vr[strongest_power_idx2];  // subcluster 2 of 2nd strongest
sns_vr[num_path + 3] = sns_vr[strongest_power_idx2];  // subcluster 3 of 2nd strongest
```

---

## 3. Per-Element 감쇠 공식: compute_sns_attenuation()

### Eq. 7.6-58 (3GPP 원문 OMML 파싱 결과)

```
α_{i,n} = 1                                    , |x_i - x_0| ≤ a AND |y_i - y_0| ≤ b  (VR 내부)
         exp(-C · |x_i - x_a| / D_n)           , |x_i - x_0| > a AND |y_i - y_0| ≤ b  (수평 이탈)
         exp(-C · |y_i - y_b| / D_n)           , |x_i - x_0| ≤ a AND |y_i - y_0| > b  (수직 이탈)
         exp(-C · sqrt(dx² + dy²) / D_n)       , otherwise                              (양방향 이탈)
```

여기서:
- `(x_i, y_i)`: 안테나 element의 위치 (wavelength 단위)
- `(x_0, y_0)`: VR reference corner
- `(x_a, y_b)`: VR 내부 경계 (reference corner 반대쪽)
- `D_n = sqrt((x_A - x_a)² + (y_B - y_b)²)`: 어레이 먼 코너까지 거리
- `C = 13`: roll-off factor (Table 7.6.14.1.2-3, 전 시나리오 동일)

**핵심**: 분수 구조는 `d / D_n`이고 `C`가 이 분수에 곱해짐 → `exp(-C · d / D_n)`.
`exp(-d / (D_n · C))`가 아님! (D_n·C를 분모에 넣으면 감쇠가 약 10배 약해짐)

### 수치 예시 (Config 3 UMi, 24×32, dH=0.5λ, dV=0.7λ)

| 위치 | d/D_n | exp(-13 · d/D_n) | dB |
|------|-------|------------------|-----|
| VR 내부 | 0 | 1.0 | 0 dB |
| VR 경계 바로 밖 | 0.05 | 0.522 | -2.8 dB |
| 어레이 중간 | 0.5 | 1.5e-3 | -28 dB |
| Far corner (d=D_n) | 1.0 | 2.3e-6 | -56 dB |

→ VR 바깥은 급격히 감쇠 (sharp visibility cutoff).

### 코드 (h/channel.h)

```cpp
static inline Real compute_sns_attenuation(
    Real pos_h, Real pos_v, const ClusterVR& vr, Real rolloff_C)
{
    if (!vr.limited) return 1.0;
    if (rolloff_C <= 0.0 || vr.D_n <= 0.0) return 1.0;

    bool inside_x = (fabs(pos_h - vr.x0) <= vr.a);
    bool inside_y = (fabs(pos_v - vr.y0) <= vr.b);
    if (inside_x && inside_y) return 1.0;

    Real dx = fabs(pos_h - vr.xa);
    Real dy = fabs(pos_v - vr.yb);

    Real d;
    if (!inside_x && inside_y)      d = dx;
    else if (inside_x && !inside_y) d = dy;
    else                            d = sqrt(dx*dx + dy*dy);

    Real gamma = exp(-rolloff_C * d / vr.D_n);  // Eq. 7.6-58
    return sqrt(gamma);  // power → amplitude
}
```

반환값은 **amplitude** (sqrt(gamma)). 채널 계수 H_usn에 직접 곱하면 됨.

---

## 4. SNS 적용 위치

### 4.1 채널 계수 생성 (GetChannelImpulseResponse)

`H_usn[n](u, s)` 계산 후 per-element 감쇠 적용:

```cpp
Real sns_atten = 1.0;
if (sns_any_limited) {
    Real pos_h = s_n * BS_dH;
    Real pos_v = s_m * BS_dV;
    sns_atten = compute_sns_attenuation(pos_h, pos_v, sns_vr[n], g_sns_rolloff_C);
}
H_usn[n](u, s) = norm * sns_atten * h_val;
```

이 부분은 per-element H_usn 레벨이므로 port-level 합산 전에 적용됨 → 올바름.

### 4.2 RSRP 계산 (find_best_tx_beam, Get_BS_antenna_field_pattern)

RSRP 공식 (TR 36.873 Eq. 8.1-1~5)에서 port-level TX 패턴은 element의 coherent sum:

```
F_{tx,p} = Σ_s w_s · exp(j·phase_s) · F_{tx,s}
```

SNS는 **coherent sum 안에서** element별로 적용해야 함:

```
F'_{tx,p} = Σ_s w_s · sqrt(α_{s,n}) · exp(j·phase_s) · F_{tx,s}
```

**주의**: 평균 감쇠를 사후에 곱하는 방식(`|F|² · avg(α)`)은 수학적으로 틀림.
`|Σ w·sqrt(α)·e^{jφ}|² ≠ |Σ w·e^{jφ}|² · avg(α)`

적용 위치:
- `find_best_tx_beam()`: 6개 AF 루프 (Handheld/Omni/Directional × NLOS/LOS)
- `Get_BS_antenna_field_pattern()`: 파라미터로 ClusterVR 포인터 전달
- `Get_antgain()` / `Get_TX_SmallScale_antgain()`: 클러스터별 VR 포인터 전달

`Get_RSRP()`의 기존 `sns_rsrp_power_atten[i]` 사후 보정은 제거됨 (per-element 방식으로 대체).

---

## 5. 글로벌 파라미터 (h/Initiallization.h)

```cpp
bool g_sns_bs_enabled;              // BS-side SNS on/off
Real g_sns_mu_P_vis;                // μ_vis (가시성 확률 평균)
Real g_sns_sigma_P_vis;             // σ_vis (가시성 확률 표준편차)
Real g_sns_vr_A, g_sns_vr_B;       // A, B (Eq. 7.6-57)
Real g_sns_vr_R;                    // R (Eq. 7.6-57 분모)
Real g_sns_vr_delta;                // σ_δ (Eq. 7.6-57, δ ~ N(0, σ_δ))
Real g_sns_rolloff_C = 13.0;       // C (Table 7.6.14.1.2-3, 전 시나리오 동일)
```

---

## 6. 자주 실수하는 포인트 (수정 이력)

### 6.1 δ 배치 (Eq. 7.6-57)

```
✗ 틀림:  V_n = A · exp(-(P_max - P_n)/R + δ) + B
✓ 올바름: V_n = A · exp(-(P_max - P_n)/R) + B + δ
```

δ는 exp 바깥에서 additive로 더해짐.

### 6.2 감쇠 공식 (Eq. 7.6-58)

```
✗ 틀림 (원래 코드):  gamma = exp(-d / (D_n · C))     → 감쇠 너무 약함 (far corner -0.33 dB)
✗ 틀림 (중간 수정):  gamma = exp(-d / C)              → D_n 정규화 누락
✓ 올바름 (OMML 파싱): gamma = exp(-C · d / D_n)       → 올바른 sharp cutoff (far corner -56 dB)
```

OMML 분수 구조: 분자 = `-d`, 분모 = `D_n`, 그리고 `·C`가 분수에 곱해지는 형태.

### 6.3 P_n 파워 (Eq. 7.6-57)

```
✗ 틀림: power[n] (K-factor 미적용, NLOS만)
✓ 올바름: powerForAngles[n] (K-factor 적용, Eq. 7.5-8)
```

LOS 조건에서 powerForAngles[0]에 K/(K+1)이 포함되어 있으므로,
첫 번째 클러스터가 자연스럽게 가장 큰 V_n (= 가장 큰 VR)을 받음.

### 6.4 LOS VR

```
✗ 틀림: LOS를 별도 클러스터로 분리하여 K/(K+1) 파워로 독립 VR 계산
✓ 올바름: sns_vr_los = sns_vr[0] (powerForAngles[0]에 이미 LOS 파워 반영)
```

### 6.5 서브클러스터 VR

```
✗ 틀림: sns_vr[num_path ~ num_path+3].limited = false (VR 없음)
✓ 올바름: 부모 클러스터의 VR을 상속 (Set_SUBCLUSTER에서 복사)
```

### 6.6 RSRP에서 SNS 적용 방식

```
✗ 틀림: avg(α) 사후 곱셈 — |F|² · avg(α_{s,n})
✓ 올바름: per-element coherent sum 안에서 적용 — |Σ w·sqrt(α)·e^{jφ}|²
```

---

## 7. 호출 순서 요약

```
Set_SmallScaleParameter()
  ├── Set_DELAY()
  ├── Set_POWER()           ← power[], powerForAngles[] 생성
  ├── Generate_VisibilityRegion()  ← sns_vr[0..num_path-1], sns_vr_los 생성
  ├── Compute_SNS_RSRP_Attenuation()  ← (legacy, 현재 미사용)
  ├── Find_Strong2Clusters()
  ├── Set_AOAAOD(), Set_ZOAZOD()
  ├── Set_InitialPhase()
  ├── Set_SUBCLUSTER()      ← sns_vr[num_path..+3] 부모 VR 복사
  └── ...

GetChannelImpulseResponse()
  └── H_usn[n](u,s) *= compute_sns_attenuation(...)  ← per-element 적용

find_best_tx_beam() / Get_BS_antenna_field_pattern()
  └── AF loop 내 sns_w = compute_sns_attenuation(...)  ← per-element 적용
```
