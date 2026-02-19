# 8 Simulations

> *Editor's note: This clause includes baseline simulation results (corresponding to a number of antenna ports and transmission scheme supported by Rel-11) with the modified evaluation methodology*

## 8.1 RSRP calculation formula

**3GPP TR 36.873 V12.7.0 (2017-12), Release 12, Page 38**

For channel model calibration and baseline performance evaluation, the following RSRP calculation formula for Tx antenna port $p$ needed for UE attachment is used (applying $K_R = 0$ for NLOS UE). Note that the TX power is CRS transmitted power per RE. The notation below is according to equation (7.3-22). $w_s$ ($s=1, \ldots, S$) represents a complex weight vector used for virtualization of port $p$ and $U$ is the number of receive antenna elements.

$$
RSRP_p = PL \cdot SF \cdot \sum_{u=1}^{U} \left( \left| \alpha_{0,u,p} \right|^2 + \sum_{n=1}^{N} \sum_{m=1}^{M} \left| \alpha_{n,m,u,p} \right|^2 \right) \cdot \frac{TX_{power}}{U}
\tag{8.1-1}
$$

where for NLOS path for $n=1, \ldots, N$,

$$
\alpha_{n,m,u,p} = \sqrt{\frac{P_n}{M(K_R+1)}}
\begin{bmatrix} F_{rx,u,\theta}(\theta_{n,m,ZOA}, \phi_{n,m,AOA}) \\ F_{rx,u,\phi}(\theta_{n,m,ZOA}, \phi_{n,m,AOA}) \end{bmatrix}^T
\begin{bmatrix} \exp(j\Phi^{\theta\theta}_{n,m}) & \sqrt{\kappa_{n,m}^{-1}} \exp(j\Phi^{\theta\phi}_{n,m}) \\ \sqrt{\kappa_{n,m}^{-1}} \exp(j\Phi^{\phi\theta}_{n,m}) & \exp(j\Phi^{\phi\phi}_{n,m}) \end{bmatrix}
\begin{bmatrix} F_{tx,p,\theta}(\theta_{n,m,ZOD}, \phi_{n,m,AOD}) \\ F_{tx,p,\phi}(\theta_{n,m,ZOD}, \phi_{n,m,AOD}) \end{bmatrix}
\tag{8.1-2}
$$

and for LOS path

$$
\alpha_{0,u,p} = \sqrt{\frac{K_R}{K_R+1}}
\begin{bmatrix} F_{rx,u,\theta}(\theta_{LOS,ZOA}, \phi_{LOS,AOA}) \\ F_{rx,u,\phi}(\theta_{LOS,ZOA}, \phi_{LOS,AOA}) \end{bmatrix}^T
\begin{bmatrix} \exp(j\Phi_{LOS}) & 0 \\ 0 & -\exp(j\Phi_{LOS}) \end{bmatrix}
\begin{bmatrix} F_{tx,p,\theta}(\theta_{LOS,ZOD}, \phi_{LOS,AOD}) \\ F_{tx,p,\phi}(\theta_{LOS,ZOD}, \phi_{LOS,AOD}) \end{bmatrix}
\tag{8.1-3}
$$

with

$$
F_{tx,p,\theta}(\theta_{n,m,ZOD}, \phi_{n,m,AOD}) = \sum_{s=1}^{S} w_s \exp\left( j 2\pi \lambda_0^{-1} \left( \hat{r}_{tx,n,m}^T \cdot \bar{d}_{tx,s} \right) \right) F_{tx,s,\theta}(\theta_{n,m,ZOD}, \phi_{n,m,AOD})
\tag{8.1-4}
$$

and

$$
F_{tx,p,\phi}(\theta_{n,m,ZOD}, \phi_{n,m,AOD}) = \sum_{s=1}^{S} w_s \exp\left( j 2\pi \lambda_0^{-1} \left( \hat{r}_{tx,n,m}^T \cdot \bar{d}_{tx,s} \right) \right) F_{tx,s,\phi}(\theta_{n,m,ZOD}, \phi_{n,m,AOD})
\tag{8.1-5}
$$

---

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
| $F_{rx,u,\theta/\phi}$ | RX 엘리먼트 $u$의 field pattern ($\theta$, $\phi$ 성분) |
| $F_{tx,s,\theta/\phi}$ | TX 엘리먼트 $s$의 element-level field pattern |
| $F_{tx,p,\theta/\phi}$ | TX port $p$의 port-level field pattern (식 8.1-4, 8.1-5) |
| $\kappa_{n,m}$ | Cross-polarization ratio (XPR) |
| $\Phi^{xx}_{n,m}$ | 랜덤 초기 위상 ($\theta\theta$, $\theta\phi$, $\phi\theta$, $\phi\phi$) |
| $\hat{r}_{tx,n,m}$ | 레이 $(n,m)$의 TX 출발 방향 단위벡터 |
| $\bar{d}_{tx,s}$ | TX 엘리먼트 $s$의 위치 벡터 |
| $\lambda_0$ | 반송파 파장 |
