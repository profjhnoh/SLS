# TR 38.901 Handheld UT 안테나 모델 상세 명세

> **Source:** 3GPP TR 38.901 V19.1.0 — Section 7.3, 7.6.14.2
> **참고 문서:** `doc/UT antenna model v2.pdf`, `doc/UT antenna model.docx`

---

## 1. 디바이스 형상 및 좌표계

### 1.1 디바이스 치수

Handheld UT는 **평면 직사각형 판**으로 모델링한다:

$$
\text{Device dimensions} = (X,\; Y,\; Z) = (15 \text{ cm},\; 7 \text{ cm},\; 0 \text{ cm}) \quad \text{(depth, width, height)}
$$

### 1.2 UT Local Coordinate System (UT LCS)

디바이스 중심이 원점이며, **기준 방향** ($\Omega_{UT}=0°$) 에서:

- $x'$축: 긴 변 방향 (15 cm, depth)
- $y'$축: 짧은 변 방향 (7 cm, width)
- $z'$축: 디바이스 평면에 수직 (법선, 화면 바깥)

기준 방향에서 GCS와 UT LCS가 일치한다 (Figure 7.3-3):

$$
\Omega_{UT,\alpha} = 0°, \quad \Omega_{UT,\beta} = 0°, \quad \Omega_{UT,\gamma} = 0°
$$

**Reference UT orientation vector**: 디바이스 평면에 수직인 $\hat{z}'$ 방향

**Near-field phase 기준점**: 디바이스 평면의 중심 (원점)

---

## 2. 안테나 후보 위치 (8개)

직사각형의 **4개 꼭짓점 + 4개 변 중앙** = 총 8개 (Figure 7.3-2).

> **Note:** CPE (9개 위치)와 달리 handheld는 중심 위치 없이 8개이다.

### 2.1 번호 체계 (Figure 7.3-2, top-down view)

Figure 7.3-3 의 bottom-left plan view ($x'$ 오른쪽, $y'$ 위쪽) 기준:

```
       y'
       ↑
  (7)---(6)---(5)       +W/2
   |           |
  (8)    ·    (4)        0    (· = 디바이스 중심)
   |           |
  (1)---(2)---(3)       -W/2
       -L/2    0   +L/2  → x'
```

번호는 좌하(1)에서 시계방향으로 순환: (1)→(2)→(3)→(4)→(5)→(6)→(7)→(8)

### 2.2 UT LCS 좌표 ($L=15$ cm, $W=7$ cm)

| 위치 | $x'$ [cm] | $y'$ [cm] | $z'$ [cm] | 유형 |
|:---:|:---:|:---:|:---:|:---:|
| **(1)** | $-L/2 = -7.5$ | $-W/2 = -3.5$ | $0$ | 꼭짓점 (좌하) |
| **(2)** | $0$ | $-W/2 = -3.5$ | $0$ | 변 중앙 (하) |
| **(3)** | $+L/2 = +7.5$ | $-W/2 = -3.5$ | $0$ | 꼭짓점 (우하) |
| **(4)** | $+L/2 = +7.5$ | $0$ | $0$ | 변 중앙 (우) |
| **(5)** | $+L/2 = +7.5$ | $+W/2 = +3.5$ | $0$ | 꼭짓점 (우상) |
| **(6)** | $0$ | $+W/2 = +3.5$ | $0$ | 변 중앙 (상) |
| **(7)** | $-L/2 = -7.5$ | $+W/2 = +3.5$ | $0$ | 꼭짓점 (좌상) |
| **(8)** | $-L/2 = -7.5$ | $0$ | $0$ | 변 중앙 (좌) |

> 구현 시 cm → m 변환: $L = 0.15$ m, $W = 0.07$ m

### 2.3 안테나 방향 벡터 (Orientation)

각 안테나는 **디바이스 중심에서 해당 안테나 위치로 향하는 방향**으로 정렬된다:

$$
\vec{v}_i = \vec{p}_i - \vec{c} = \vec{p}_i, \qquad
\hat{u}_i = \frac{\vec{v}_i}{\|\vec{v}_i\|}
$$

| 위치 | 중심→안테나 $\vec{v}_i$ | 정규화 $\hat{u}_i$ | 방위각 $\phi_c$ |
|:---:|:---:|:---:|:---:|
| (1) | $(-7.5,\; -3.5,\; 0)$ | $(-0.906,\; -0.423,\; 0)$ | $\approx -155°$ |
| (2) | $(0,\; -3.5,\; 0)$ | $(0,\; -1,\; 0)$ | $-90°$ |
| (3) | $(+7.5,\; -3.5,\; 0)$ | $(+0.906,\; -0.423,\; 0)$ | $\approx -25°$ |
| (4) | $(+7.5,\; 0,\; 0)$ | $(+1,\; 0,\; 0)$ | $0°$ |
| (5) | $(+7.5,\; +3.5,\; 0)$ | $(+0.906,\; +0.423,\; 0)$ | $\approx +25°$ |
| (6) | $(0,\; +3.5,\; 0)$ | $(0,\; +1,\; 0)$ | $+90°$ |
| (7) | $(-7.5,\; +3.5,\; 0)$ | $(-0.906,\; +0.423,\; 0)$ | $\approx +155°$ |
| (8) | $(-7.5,\; 0,\; 0)$ | $(-1,\; 0,\; 0)$ | $180°$ |

---

## 3. 단일 엘리먼트 방사 패턴 (Table 7.3-2)

### Table 7.3-2: Radiation power pattern of a single antenna element for handheld UT

안테나는 $\hat{\theta}''$, $\hat{\varphi}''$ 방향으로 oriented 되어 있으며, boresight는 $\theta'' = 90°$, $\varphi'' = 0°$ 방향 ($x'$축, broadside)이다.

#### 수직 컷 (Vertical cut)

$$
A''_{E,V}(\theta'') = -\min\!\left\{ 12 \left(\frac{\theta'' - 90°}{\theta_{3dB}}\right)^2,\; SLA_V \right\} \quad \text{[dB]}
$$

$$
\theta_{3dB} = 125°, \quad SLA_V = 22.5 \text{ dB}, \quad \theta'' \in [0°,\; 180°]
$$

#### 수평 컷 (Horizontal cut)

$$
A''_{E,H}(\varphi'') = -\min\!\left\{ 12 \left(\frac{\varphi''}{\varphi_{3dB}}\right)^2,\; A_{\max} \right\} \quad \text{[dB]}
$$

$$
\varphi_{3dB} = 125°, \quad A_{\max} = 22.5 \text{ dB}, \quad \varphi'' \in [-180°,\; 180°]
$$

#### 3D 방사 패턴

$$
A''_{E}(\theta'', \varphi'') = -\min\!\left\{ -\Big(A''_{E,V}(\theta'') + A''_{E,H}(\varphi'')\Big),\; A_{\max} \right\} \quad \text{[dB]}
$$

#### 최대 지향성 이득

$$
G_{E,\max} = 5.3 \text{ dBi}
$$

#### 최종 엘리먼트 패턴

$$
A'(\theta'', \varphi'') = G_{E,\max} + A''_E(\theta'', \varphi'') \quad \text{[dBi]}
$$

linear scale 변환:

$$
A'_{\text{linear}}(\theta'', \varphi'') = 10^{A'(\theta'',\varphi'')/10}
$$

#### NOTE

> For UT antenna modelling of handheld devices, optional antenna imbalance can be modelled. If modelled, randomized loss is applied per UT antenna port and randomized loss can be applied independently for the UL and DL directions. No imbalance is modelled by default.

### Table 7.3-1 (BS)과의 비교

| 파라미터 | Table 7.3-2 (Handheld UT) | Table 7.3-1 (BS) |
|:---:|:---:|:---:|
| $\theta_{3dB}$ | $125°$ | $65°$ |
| $\varphi_{3dB}$ | $125°$ | $65°$ |
| $SLA_V$ | $22.5$ dB | $30$ dB |
| $A_{\max}$ | $22.5$ dB | $30$ dB |
| $G_{E,\max}$ | $5.3$ dBi | $8$ dBi |

Handheld UT는 BS 대비 빔폭이 약 2배 넓고, 최대 감쇄가 낮으며, 이득이 작다.

---

## 4. 편파 안테나 모델링 (Section 7.3.2: Polarized Antenna Modelling)

> 이 섹션은 TR 38.901 Section 7.3.2의 **Handheld UT Model** 부분을 상세하게 정리한다.

### 4.1 Handheld UT 모델 — 기준 방사 패턴 (Reference Radiation Pattern)

기준 좌표계(LCS'')에서 안테나 방사 패턴의 편파 성분을 정의한다. 기준 편파 방향은 $\hat{z}''$ 축 방향 ($\hat{\theta}''$ 편파, 수직 편파)이다.

#### Case 1: Single Polarization (1개 안테나 패턴 — 주로 FR1)

하나의 후보 안테나 위치에 **하나의 안테나 field pattern**을 사용하는 경우:

$$
F''_{\theta''}(\theta'', \varphi'') = \sqrt{A'_{\text{linear}}(\theta'', \varphi'')}, \qquad F''_{\varphi''}(\theta'', \varphi'') = 0
$$

기준 패턴은 $\hat{\theta}''$ 방향으로만 이득을 가지며 ($\hat{z}''$ 축 편파), 이를 **"polarization direction along the $Z''$ axis"** 라 한다.

#### Case 2: Dual Polarization (2개 안테나 패턴 — FR2 전용, not intended for FR1)

하나의 후보 안테나 위치에 **두 개의 서로 다른 편파 field pattern**을 사용하는 경우:

$$
\text{첫 번째 패턴: } \quad F''_{\theta''}(\theta'', \varphi'') = \sqrt{A'_{\text{linear}}(\theta'', \varphi'')}, \quad F''_{\varphi''}(\theta'', \varphi'') = 0
$$

$$
\text{두 번째 패턴: } \quad F''_{\theta''}(\theta'', \varphi'') = 0, \quad F''_{\varphi''}(\theta'', \varphi'') = \sqrt{A'_{\text{linear}}(\theta'', \varphi'')}
$$

### 4.2 Stage 1: 안테나 개별 회전 (Eq. 7.3-6, 7.3-7, 7.3-8)

각 편파 field component를 각 UT 안테나 $u$의 **방향(orientation)**과 **편파 방향(polarization direction)**에 따라 회전하여, UT LCS 좌표계의 field pattern $(F'_{u,\theta'}, F'_{u,\varphi'})$를 구한다.

#### Eq. 7.3-6: 편파 회전 행렬

$$
\begin{bmatrix}
F'_{u,\theta'}(\theta', \varphi') \\
F'_{u,\varphi'}(\theta', \varphi')
\end{bmatrix}
=
\begin{bmatrix}
+\cos\psi_u & -\sin\psi_u \\
+\sin\psi_u & +\cos\psi_u
\end{bmatrix}
\begin{bmatrix}
F''_{\theta''}(\theta'', \varphi'') \\
F''_{\varphi''}(\theta'', \varphi'')
\end{bmatrix}
\tag{7.3-6}
$$

여기서 $(\theta', \varphi')$는 UT LCS 좌표, $(\theta'', \varphi'')$는 안테나 개별 LCS'' 좌표이다.

#### Eq. 7.3-7: $\cos\psi_u$ 정의

$$
\cos\psi_u = \frac{
\cos\beta_u \cos\gamma_u \sin\theta'
- \Big(\sin\beta_u \cos\gamma_u \cos(\varphi'-\alpha_u) - \sin\gamma_u \sin(\varphi'-\alpha_u)\Big) \cos\theta'
}{
\sqrt{1 - \Big(\cos\beta_u \cos\gamma_u \cos\theta' + \big(\sin\beta_u \cos\gamma_u \cos(\varphi'-\alpha_u) - \sin\gamma_u \sin(\varphi'-\alpha_u)\big) \sin\theta'\Big)^2}
}
\tag{7.3-7}
$$

#### Eq. 7.3-8: $\sin\psi_u$ 정의

$$
\sin\psi_u = \frac{
\sin\beta_u \cos\gamma_u \sin(\varphi'-\alpha_u) + \sin\gamma_u \cos(\varphi'-\alpha_u)
}{
\sqrt{1 - \Big(\cos\beta_u \cos\gamma_u \cos\theta' + \big(\sin\beta_u \cos\gamma_u \cos(\varphi'-\alpha_u) - \sin\gamma_u \sin(\varphi'-\alpha_u)\big) \sin\theta'\Big)^2}
}
\tag{7.3-8}
$$

> **핵심:** Eq. 7.3-7, 7.3-8에서의 입력 각도는 **UT LCS 좌표 $(\theta', \varphi')$** 이다. 이는 이전 섹션의 Eq. 7.1-15 (GCS 입력)와 **다른 좌표계**에서 동작한다.

#### 각도 매핑 (Eq. 7.1-7, 7.1-8 적용)

UT LCS 각도 $(\theta', \varphi')$를 안테나 $u$의 LCS'' 각도 $(\theta'', \varphi'')$로 변환한다. 이때 3D 회전 각도 $(\alpha_u, \beta_u, \gamma_u)$를 사용한다:

$$
\theta'' = \arccos\!\Big(
\cos\beta_u \cos\theta'
+ \sin\beta_u \sin\theta' \cos(\varphi' - \alpha_u) \cos\gamma_u
- \sin\beta_u \sin\theta' \sin(\varphi' - \alpha_u) \sin\gamma_u
\Big)
\tag{7.1-7 적용}
$$

$$
\varphi'' = \arg\!\Big(
\big[\cos\beta_u \sin\theta' \cos(\varphi'-\alpha_u)\cos\gamma_u + \sin\beta_u\cos\theta'\cos\gamma_u
- \cos\beta_u \sin\theta' \sin(\varphi'-\alpha_u)\sin\gamma_u\big]
+ j\big[\cos\gamma_u \sin\theta'\sin(\varphi'-\alpha_u) + \sin\gamma_u\cos\theta'\big]
\Big)
\tag{7.1-8 적용}
$$

> **Note:** 여기서 Eq. 7.1-7, 7.1-8의 $(\alpha, \beta, \gamma)$에 안테나별 $(\alpha_u, \beta_u, \gamma_u)$를 대입하고, 입력/출력 각도에 UT LCS ↔ LCS''를 사용한다.

#### $\gamma_u = 0$ 일 때 단순화 (Single Polarization Handheld)

$$
\cos\psi_u = \frac{
\cos\beta_u \sin\theta' - \sin\beta_u \cos(\varphi'-\alpha_u) \cos\theta'
}{
\sqrt{1 - \Big(\cos\beta_u \cos\theta' + \sin\beta_u \cos(\varphi'-\alpha_u) \sin\theta'\Big)^2}
}
$$

$$
\sin\psi_u = \frac{
\sin\beta_u \sin(\varphi'-\alpha_u)
}{
\sqrt{1 - \Big(\cos\beta_u \cos\theta' + \sin\beta_u \cos(\varphi'-\alpha_u) \sin\theta'\Big)^2}
}
$$

$\beta_u = 90°$, $\gamma_u = 0°$ 추가 적용 시:

$$
\cos\psi_u = \frac{
-\sin(\varphi'-\alpha_u) \cos\theta'  \cdot \sin 90°
}{
\sqrt{1 - \Big(\cos(\varphi'-\alpha_u) \sin\theta'\Big)^2}
}
= \frac{-\sin(\varphi'-\alpha_u)\cos\theta'}{\sqrt{1-\cos^2(\varphi'-\alpha_u)\sin^2\theta'}}
$$

$$
\sin\psi_u = \frac{
\sin(\varphi'-\alpha_u)  \cdot \sin 90°
}{
\sqrt{1 - \Big(\cos(\varphi'-\alpha_u) \sin\theta'\Big)^2}
}
= \frac{\sin(\varphi'-\alpha_u)}{\sqrt{1-\cos^2(\varphi'-\alpha_u)\sin^2\theta'}}
$$

### 4.3 안테나 위치별 편파 방향 — Single Polarization (Figure 7.3-7)

각 안테나의 실제 편파 방향은 다음 두 조건을 만족한다:

1. 디바이스 **평면에 평행** (Parallel with the plane of the handheld UT)
2. 중심→안테나 방향 벡터에 **수직** (Perpendicular to the direction from UT center to candidate antenna location)

중심→안테나 방향 벡터 $\vec{d}_u = (dx_u,\; dy_u)$ 에 대해, single-pol 편파 방향:

$$
\vec{p}_u = \frac{(-dy_u,\; dx_u,\; 0)}{\sqrt{dx_u^2 + dy_u^2}}
$$

| 위치 | 중심→안테나 $(dx, dy)$ | 편파 방향 $\vec{p}_u$ | 편파 방위각 |
|:---:|:---:|:---:|:---:|
| (1) | $(-7.5,\; -3.5)$ | $(+3.5,\; -7.5,\; 0)/8.28$ | $\approx -65°$ |
| (2) | $(0,\; -3.5)$ | $(+1,\; 0,\; 0)$ | $0°$ |
| (3) | $(+7.5,\; -3.5)$ | $(+3.5,\; +7.5,\; 0)/8.28$ | $\approx +65°$ |
| (4) | $(+7.5,\; 0)$ | $(0,\; +1,\; 0)$ | $+90°$ |
| (5) | $(+7.5,\; +3.5)$ | $(-3.5,\; +7.5,\; 0)/8.28$ | $\approx +115°$ |
| (6) | $(0,\; +3.5)$ | $(-1,\; 0,\; 0)$ | $180°$ |
| (7) | $(-7.5,\; +3.5)$ | $(-3.5,\; -7.5,\; 0)/8.28$ | $\approx -115°$ |
| (8) | $(-7.5,\; 0)$ | $(0,\; -1,\; 0)$ | $-90°$ |

> Figure 7.3-7에서 각 안테나 위치의 화살표가 이 편파 방향을 나타낸다.

### 4.4 안테나 위치별 편파 방향 — Dual Polarization (Figure 7.3-8, FR2 전용)

하나의 후보 안테나 위치에 두 개의 서로 다른 편파 field pattern을 사용하는 경우 (not intended for FR1):

- **첫 번째 패턴**: Figure 7.3-7의 single-pol 편파 방향에서, **중심→안테나 방향 벡터를 회전축으로 $+45°$ 추가 회전** (rotated 45 degrees about the direction from UT center to candidate antenna location)
- **두 번째 패턴**: 첫 번째 패턴의 편파 방향에 **수직**이고, 중심→안테나 방향에도 **수직** (perpendicular to the polarization direction of the first antenna field pattern and perpendicular to the direction from UT center to candidate antenna location)

> Figure 7.3-8은 위치 (6)을 예시로 한 측면도(side view)에서 dual-pol 편파 방향을 도시한다. 중심→안테나 방향 벡터를 축으로, 기준 편파를 ±45° 기울인 형태이다.

### 4.5 3D 회전 각도 $(\alpha_u, \beta_u, \gamma_u)$의 결정

Eq. 7.3-6~8의 $(\alpha_u, \beta_u, \gamma_u)$는 **각 UT 안테나 $u$의 방향(orientation)과 편파 방향(polarization direction)**에 따라 결정된다.

기준 패턴의 $\hat{\theta}''$ 편파(수직, $\hat{z}''$ 방향)를 각 안테나의 실제 편파 방향으로 변환하기 위한 3D 회전:

$$
\alpha_u = \text{atan2}(p_{u,y},\; p_{u,x}) \quad \text{(편파 방향의 방위각)}
$$

$$
\beta_u = 90° \quad \text{(}\hat{z}''\text{ 방향 편파를 디바이스 평면으로 눕힘)}
$$

$$
\gamma_u = 0° \quad \text{(single polarization — 추가 슬랜트 없음)}
$$

| 위치 | $\alpha_u$ | $\beta_u$ | $\gamma_u$ |
|:---:|:---:|:---:|:---:|
| (1) | $\text{atan2}(-7.5,\; 3.5) \approx -65°$ | $90°$ | $0°$ |
| (2) | $\text{atan2}(0,\; 1) = 0°$ | $90°$ | $0°$ |
| (3) | $\text{atan2}(7.5,\; 3.5) \approx 65°$ | $90°$ | $0°$ |
| (4) | $\text{atan2}(1,\; 0) = 90°$ | $90°$ | $0°$ |
| (5) | $\text{atan2}(7.5,\; -3.5) \approx 115°$ | $90°$ | $0°$ |
| (6) | $\text{atan2}(0,\; -1) = 180°$ | $90°$ | $0°$ |
| (7) | $\text{atan2}(-7.5,\; -3.5) \approx -115°$ | $90°$ | $0°$ |
| (8) | $\text{atan2}(-1,\; 0) = -90°$ | $90°$ | $0°$ |

> Dual polarization의 경우 $\gamma_u = \pm 45°$ 로 변경하여 추가 슬랜트를 적용한다.

### 4.6 Stage 2: UT 전체 GCS 회전 (Eq. 7.1-11)

Stage 1에서 구한 UT LCS field pattern $(F'_{u,\theta'}, F'_{u,\varphi'})$를 UT의 전체 방향 $(\Omega_{UT,\alpha}, \Omega_{UT,\beta}, \Omega_{UT,\gamma})$에 따라 GCS로 최종 회전한다. 이는 Clause 7.1.3의 Eq. 7.1-11을 사용한다:

$$
\begin{bmatrix}
F_{u,\theta}(\theta, \varphi) \\
F_{u,\varphi}(\theta, \varphi)
\end{bmatrix}
=
\begin{bmatrix}
\cos\psi_{UT} & -\sin\psi_{UT} \\
\sin\psi_{UT} & \cos\psi_{UT}
\end{bmatrix}
\begin{bmatrix}
F'_{u,\theta'}(\theta', \varphi') \\
F'_{u,\varphi'}(\theta', \varphi')
\end{bmatrix}
\tag{7.1-11}
$$

$\psi_{UT}$는 $(\Omega_{UT,\alpha}, \Omega_{UT,\beta}, \Omega_{UT,\gamma})$와 GCS 각도 $(\theta, \varphi)$로 Eq. 7.1-15를 적용하여 구한다:

$$
\psi_{UT} = \text{atan2}\!\left(
\sin\Omega_\beta \sin\Omega_\gamma \cos(\varphi-\Omega_\alpha) + \cos\Omega_\beta \sin\Omega_\gamma \sin\theta - \sin\Omega_\beta \cos\Omega_\gamma \cos\theta,\;
\cos\Omega_\beta \cos\Omega_\gamma - \sin\Omega_\beta \sin\Omega_\gamma \cos(\varphi-\Omega_\alpha) \sin\theta
\right)
\tag{7.1-15}
$$

> 약칭: $\Omega_\alpha = \Omega_{UT,\alpha}$, $\Omega_\beta = \Omega_{UT,\beta}$, $\Omega_\gamma = \Omega_{UT,\gamma}$

UT LCS 각도 $(\theta', \varphi')$와 GCS 각도 $(\theta, \varphi)$ 간의 매핑에도 동일하게 Eq. 7.1-7, 7.1-8을 $(\Omega_{UT,\alpha}, \Omega_{UT,\beta}, \Omega_{UT,\gamma})$로 적용한다.

### 4.7 전체 변환 요약 (3개 좌표계)

$$
\underbrace{(F''_{\theta''},\; F''_{\varphi''})}_{\text{기준 패턴 (LCS'')}}
\;\xrightarrow[\text{Eq. 7.3-6/7/8}]{\text{Stage 1: } (\alpha_u, \beta_u, \gamma_u)}\;
\underbrace{(F'_{u,\theta'},\; F'_{u,\varphi'})}_{\text{UT LCS}}
\;\xrightarrow[\text{Eq. 7.1-11}]{\text{Stage 2: } (\Omega_{UT,\alpha}, \Omega_{UT,\beta}, \Omega_{UT,\gamma})}\;
\underbrace{(F_{u,\theta},\; F_{u,\varphi})}_{\text{GCS}}
$$

| 단계 | 좌표계 변환 | 회전 파라미터 | 수식 |
|:---:|:---:|:---:|:---:|
| Stage 1 | LCS'' → UT LCS | $(\alpha_u, \beta_u, \gamma_u)$ per antenna | Eq. 7.3-6, 7.3-7, 7.3-8, 7.1-7, 7.1-8 |
| Stage 2 | UT LCS → GCS | $(\Omega_{UT,\alpha}, \Omega_{UT,\beta}, \Omega_{UT,\gamma})$ per UT | Eq. 7.1-11, 7.1-7, 7.1-8, 7.1-15 |

> **현재 코드와의 차이:** `Get_UE_antenna_pattern()`은 Stage 2 (UT→GCS) 회전만 수행한다. Handheld 모델에서는 **Stage 1 (안테나별 개별 회전)**이 추가되어 **총 2단계 회전**이 필요하다.

> **주의:** Stage 1의 $\psi_u$ 계산은 Eq. 7.3-7, 7.3-8을 사용하며, 입력 각도가 **UT LCS $(\theta', \varphi')$** 이다. Stage 2의 $\psi_{UT}$ 계산은 Eq. 7.1-15를 사용하며, 입력 각도가 **GCS $(\theta, \varphi)$** 이다. 두 수식의 형태는 유사하지만 좌표계가 다르므로 혼동하지 않도록 한다.

---

## 5. 안테나 위치 벡터와 채널 계수

### 5.1 위치 벡터 $\vec{r}_u$

채널 계수 (Eq. 7.5-22)의 phase term에서 사용:

$$
\exp\!\left(j \frac{2\pi}{\lambda_0} \hat{r}_{rx,n,m} \cdot \vec{r}_u\right)
$$

각 안테나의 위치 벡터 $\vec{r}_u$ 는 디바이스 중심 기준 UT LCS 좌표이며, UT의 GCS 회전에 따라 변환:

$$
\vec{r}_u^{GCS} = R_{UT}(\Omega_{UT,\alpha}, \Omega_{UT,\beta}, \Omega_{UT,\gamma}) \cdot \vec{r}_u^{LCS}
$$

| 위치 | $\vec{r}_u^{LCS}$ [m] |
|:---:|:---:|
| (1) | $(-0.075,\; -0.035,\; 0)$ |
| (2) | $(0,\; -0.035,\; 0)$ |
| (3) | $(+0.075,\; -0.035,\; 0)$ |
| (4) | $(+0.075,\; 0,\; 0)$ |
| (5) | $(+0.075,\; +0.035,\; 0)$ |
| (6) | $(0,\; +0.035,\; 0)$ |
| (7) | $(-0.075,\; +0.035,\; 0)$ |
| (8) | $(-0.075,\; 0,\; 0)$ |

> 현재 코드의 `d_rx[M][N][P][Mg][Ng]` 배열에 해당. Handheld 모델에서는 균일 배열 대신 **8개 후보 위치의 실제 좌표**를 사용해야 한다.

### 5.2 UT 회전 행렬

$$
R_{UT} = R_z(\Omega_{UT,\alpha}) \cdot R_y(\Omega_{UT,\beta}) \cdot R_x(\Omega_{UT,\gamma})
$$

$$
R_x(\gamma) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & \cos\gamma & -\sin\gamma \\ 0 & \sin\gamma & \cos\gamma \end{bmatrix}, \quad
R_y(\beta) = \begin{bmatrix} \cos\beta & 0 & \sin\beta \\ 0 & 1 & 0 \\ -\sin\beta & 0 & \cos\beta \end{bmatrix}, \quad
R_z(\alpha) = \begin{bmatrix} \cos\alpha & -\sin\alpha & 0 \\ \sin\alpha & \cos\alpha & 0 \\ 0 & 0 & 1 \end{bmatrix}
$$

---

## 6. Calibration 방향

| 시나리오 | $\Omega_{UT,\alpha}$ | $\Omega_{UT,\beta}$ | $\Omega_{UT,\gamma}$ |
|:---:|:---:|:---:|:---:|
| One-hand blockage (SNS) | $\text{uniform}[0°, 360°]$ | $45°$ | $0°$ |
| Dual-hand blockage (SNS) | $\text{uniform}[0°, 360°]$ | $0°$ | $45°$ |
| Hand + head blockage (SNS) | $\text{uniform}[0°, 360°]$ | $90°$ | $0°$ |
| All other cases (default) | $\text{uniform}[0°, 360°]$ | $45°$ | $0°$ |

---

## 7. Spatial Non-Stationarity (SNS) — Section 7.6.14.2

### 7.1 적용 비율

UT의 **90%**에 SNS를 적용한다 (10%는 균일 안테나 동작 가정).

### 7.2 Usage Scenario 확률

각 UT에 대해 Table 7.6.14.2-1의 확률 분포에 따라 usage scenario를 결정:
- One-hand blockage
- Dual-hand blockage
- Hand + head blockage

### 7.3 위치별 Power Attenuation

Table 7.6.14.2-2에서 각 후보 안테나 위치 (1)~(8)에 대한 power attenuation factor $\Gamma_{RX,u}$ (수신 시) 또는 $\Gamma_{TX,u}$ (송신 시)를 결정.

### 7.4 채널 계수 적용

Step 11 (Eq. 7.5-28, 7.5-29)에서 amplitude에 $\sqrt{\Gamma}$를 곱한다:

$$
H'_{u,s,\text{NLOS}} = \sqrt{\Gamma_{RX,u}} \cdot H_{u,s,\text{NLOS}}
$$

$$
H'_{u,s,\text{LOS}} = \sqrt{\Gamma_{RX,u}} \cdot H_{u,s,\text{LOS}}
$$

> $\Gamma$는 power 비율이므로, field amplitude에는 $\sqrt{\Gamma}$를 곱한다.

---

## 8. 구현 체크리스트

### Phase 1: Element Pattern (Table 7.3-2)

- [ ] `ue_antenna_element_gain` 값을 $5.3$ dBi로 설정
- [ ] `Get_UE_antenna_pattern()` 내 $\theta_{3dB}$, $\varphi_{3dB}$ 를 $125°$로 변경
- [ ] $SLA_V$, $A_{\max}$를 $22.5$ dB로 변경

### Phase 2: 안테나 배치 및 편파

- [ ] MS 구조체에 8개 후보 위치 좌표 추가
- [ ] 안테나별 $(\alpha_u, \beta_u, \gamma_u)$ 계산 및 저장
- [ ] 안테나별 편파 방향 벡터 $\vec{p}_u$ 계산

### Phase 3: 2단계 회전

- [ ] Stage 1 회전 구현 (Eq. 7.3-6~8): 안테나 개별 $(\alpha_u, \beta_u, \gamma_u)$
- [ ] Stage 2 회전 구현 (Eq. 7.1-11): UT 전체 $(\Omega_{UT,\alpha}, \Omega_{UT,\beta}, \Omega_{UT,\gamma})$
- [ ] 기존 `Get_UE_antenna_pattern()` 에서 handheld 분기 추가

### Phase 4: 채널 계수

- [ ] `d_rx` 배열에 8개 후보 위치의 실제 좌표 설정
- [ ] 위치 벡터의 GCS 회전 적용

### Phase 5: SNS (선택)

- [ ] 90% UT에 SNS 적용 로직 추가
- [ ] Usage scenario 확률 테이블 구현
- [ ] 위치별 $\Gamma_{RX,u}$ 테이블 구현
- [ ] 채널 계수에 $\sqrt{\Gamma}$ 적용

### Phase 6: Config 파일

- [ ] Handheld 모델 on/off 파라미터 추가
- [ ] 사용할 안테나 위치 선택 파라미터 추가
