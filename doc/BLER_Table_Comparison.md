# BLER Table 비교 — 우리 BLER_5G.dat vs MATLAB 5G Toolbox L2SM

> **⚠️ 복구 노트 (2026-07-10)**: 이 문서의 원본(git 미추적 상태)이 tmp 파일 일괄 정리에
> 휩쓸려 삭제되어 세션 기록에서 복구함. 2026-07-03 이후의 모든 섹션(구현 완료 이하)은
> 트랜스크립트의 append 원문 그대로이며, 그 이전의 §1~4(원본 테이블 비교 분석)는 아래에
> **축약 재구성**함. 원본 전문이 필요하면 Dropbox 웹 버전 기록에서
> doc/BLER_Table_Comparison.md 복원 가능.

> **목적**: 시뮬레이터가 사용하는 BLER 테이블(`BLER_5G.dat`)의 정확성을 MATLAB R2026a
> 5G Toolbox의 표준 BLER 테이블(`L2SM.mat`)과 비교 검증.

## 1~4. 원본 비교 분석 [축약 재구성]

- **우리 테이블**: `BLER_5G.dat` = SNR(−9..19.5dB @0.02 step, 1426행) × CQI(15) 단일
  2-D 테이블, 값=BLER 0~1. 조회는 `bler_lookup_dB` (floor, 무보간). 코드레이트·TBS 축 없음.
- **MATLAB 테이블**: `toolbox/5g/netmod/+nr5g/+internal/L2SM.mat`의 `awgnTables` —
  (BGN, R×1024, Qm, Zc)별 [유효SINR_dB, 코드블록 BLER] 커브, BGN1 3표+BGN2 2표
  총 ~13,000 커브(Qm≤8 필터 후 17,712 export). TBS는 nrDLSCHInfo 세그먼테이션(TS
  38.212 §5.2.2 → BGN/Zc/C)을 거쳐 반영, TB BLER = 1−(1−cb)^C, log10(BLER) 선형 보간.
- **비교 결과**: 커브 형상은 유사하나 MATLAB 커브가 훨씬 가파름(~1dB LDPC 워터폴).
  같은 (변조,레이트)에서 TBS에 따라 BLER이 크게 변동(오라클: MCS4@4dB, TBS 292 vs
  3000 → ~1700×). 우리 테이블의 BLER<0.1 임계 SNR은 NR 권장 동작점과 CQI별 수 dB
  이내 차이.
- **원본 결론(당시)**: 시스템 레벨 추상화로는 우리 테이블도 충분 — 단 (코드레이트,
  TBS) 의존성 반영이 필요하면 MATLAB 테이블 도입(후일 기관 향후계획 항목으로 실행됨).


---

# 구현 완료 (2026-07-03): MATLAB L2SM 테이블 통합

위 비교의 결론("현행 유지")을 뒤집고 MATLAB 테이블을 실제 도입했다. 기관 향후계획 항목 이행.

## 아키텍처
- **Export**: `export_matlab_bler.m` → `matlab_l2sm_bler.dat` (17,712 커브, 6.8MB; Qm∈{2,4,6,8}, (BGN,R×1024,Qm,Zc)별 [SNR,cbBLER]) + 오라클 `bler_compare/matlab_reference_tuples.csv` (448 튜플, 실제 MATLAB interpolatePER/nrDLSCHInfo로 생성).
- **C++ 모듈**: `src/nr_l2sm.cpp` + `h/nr_l2sm.h` — TS 38.212 §5.2.2 세그먼테이션(nrDLSCHInfo 동등), getAWGNTable 재현(버킷→R상향→Qm/Zc정확), interpolatePER 재현(첫 0행 절단, log10 선형보간±외삽, <1e-6→0), `TB_BLER=1-(1-cb)^C`. `nr_mcs_to_qm_r`로 CR_X1024 3중복 통합.
- **플래그** (분리 A/B 가능): `matlab_bler`(BLER 코인), `matlab_cqi_thresholds`(SINR_threshold_dB 재생성), `matlab_bler_selftest`(오라클 검증). 모두 기본 0 — off 시 레거시 경로 그대로.
- **호출부 4곳 전환**: BLER_Mapping per-layer/TB, AddThroughput IR, Receive_DL 인라인(→헬퍼; 단 Receive_DL은 호출부 없는 죽은 코드로 확인됨).
- **임계값 재생성**: CQI별 매핑된 MCS의 (Qm,R) + 전대역 1레이어 기준 TBS로 TB-BLER=0.1 SNR 이분탐색. 결과: 기존 대비 +3~8dB (매핑 MCS가 CQI 명목 레이트보다 공격적 + TB단위 기준이라 자기일관적 보수화). 예: CQI4 −1.90→+3.34dB, CQI15 +19.24→+25.86dB.

## 검증
- **셀프테스트 448/448 통과** (세그먼테이션 bgn/Zc/C 정확일치, BLER log10 오차 <1e-3, 0/1 클램프 정확).
- **오라클이 TBS 의존성 실증**: MCS4@4dB에서 TBS 292 → BLER 2.0e-2 vs TBS 3000 → 1.2e-5 (~1700×). C>1 법칙 확인 (C=24: cb 0.0075→tb 0.166).
- **회귀**: 같은-바이너리 단일스레드 재현성 완벽(bit-identical). 단, **바이너리 간 bit-identity는 원리적으로 불가** — 참조 코드에 전역만 있는 더미 TU를 추가해도 4개 UE 궤적이 변함(레이아웃 민감 읽기 존재 증명; 미초기화/전역OOB 계열). flag-off 로직의 값-동일성은 코드검증(모든 CR 상수 float 정확표현)으로 별도 입증.

## 리스크 요약 (계획서 R1)
유효 SINR 정의 불일치(우리 EESM vs 커브 축 RBIR): OLLA가 평균 편향 흡수 + 임계값이 같은 EESM-조회로 도출되어 자기일관. 필요 시 Beta 재캘리브레이션 → RBIR 포팅 순 에스컬레이션.

## ⚠️ 라이선스
`matlab_l2sm_bler.dat`·오라클 CSV는 MathWorks 소유 데이터 — **연구 repo 전용, EMIMO_SLS_dist 배포 금지** (export 스크립트만 배포).

*(A/B 결과는 실행 완료 후 아래에 추가)*


## A/B 결과 (2026-07-04, UMi 8-port mx12 UE20, 2drop×500slot)

| arm | SE | eSINR | 1st-tx BLER | avgRI | p05 edge |
|-----|----:|----:|----:|----:|----:|
| legacy (matlab_bler=0) | 6.47 | 5.7dB | 0.123 | 3.88 | 0.589 |
| matlab_bler=1 | 2.33 | 4.5dB | 0.228 | 3.91 | 0.022 |
| matlab_bler=1 + thresholds | 2.46 | 4.5dB | 0.198 | 3.82 | 0.003 |

**해석 — R1 리스크(유효 SINR 정의 불일치) 발현:**
- MATLAB AWGN 커브는 LDPC 워터폴이 가파름(~1dB) → 커브 x축(RBIR-유효SINR)과 우리 EESM(레거시 Beta) 사이의 ~1-2dB 편향이 BLER을 크게 튀김.
- 레거시 Beta[15]는 구 테이블에 맞춰진 값 — AWGN-커브용 EESM 캘리브레이션이 아님.
- OLLA 단일 오프셋(min −10dB 클램프)로는 레이어/UE별 편향을 못 잡음 → BLER 0.2 수준 고착, 링크어댑테이션 하향 → SE −64%.
- 임계값 재생성(arm3)은 방향은 옳으나(BLER 0.228→0.198) OLLA 포화로 효과 제한.

**결론**: 통합 자체는 정확(셀프테스트 448/448 = 테이블 충실 재현). 새 테이블을 기본값으로 쓰려면 계획서의 보정 사다리 다음 단계가 필요:
1. **Beta[] per-MCS 재캘리브레이션** (MATLAB AWGN 커브 기준 오프라인 피팅) ← 다음 작업
2. 그래도 부족하면 RBIR ESM 포팅 (rbir.mat)

**현 상태 권장**: `matlab_bler=0` (기본값) 유지. 플래그는 연구/보정용으로 사용.


## RBIR ESM 포팅 (2026-07-04)

R1 보정 사다리 2단계 실행: 유효 SINR을 커브의 원래 x축 정의(RBIR)로 계산.

- **Export**: `export_rbir_esm.m` → `rbir_esm.dat` (mod 4/16/64/256, [SNR_dB,RBIR] + min/max 클램프) + 오라클 32튜플(실제 `wireless.internal.L2SM.calculateEffectiveSINR`로 생성).
- **알고리즘 (dB 도메인, 검증)**: 각 SINR을 [minSNR,maxSNR] 클램프 → RBIR 순매핑(선형보간) → 평균 → 역매핑. Chase 결합이 dB를 더하는 것으로 단위 확정.
- **C++**: nr_l2sm에 `rbir_effective_sinr_linear/dB` + `matlab_rbir` 플래그(기본 1, matlab_bler=1일 때 per-layer·fallback ESM을 EESM→RBIR 전환). **셀프테스트 32/32 통과** (<1e-3 dB).

### A/B 재실행 (RBIR, 동일 조건)

| arm | SE | 1st-tx BLER | (EESM이던 직전) |
|-----|----:|----:|----:|
| legacy | 5.90 | 0.134 | 6.47 / 0.123 |
| matlab_bler (RBIR) | **3.51** | 0.224 | 2.33 / 0.228 |
| + thresholds | 2.55 | 0.180 | 2.46 / 0.198 |

**해석**:
- RBIR로 matlab arm SE **+51%** (2.33→3.51) — ESM 축 불일치가 실제 편향원이었음을 확인.
- 그러나 BLER ~0.22로 타깃 0.1 대비 갭 잔존. 남은 후보: ① 500슬롯 중 OLLA 적응 400슬롯뿐(전이 미수렴; 1200슬롯 진단으로 검증), ② per-RB 입도(MATLAB은 per-subcarrier), ③ AWGN-가파른 커브가 per-RB SINR 모델 오차를 직격 — 레거시 테이블은 자체 생태계(Beta·임계값·OLLA)에 end-to-end 튜닝되어 있었음.
- **임계값 재생성이 RBIR 하에선 역효과**(2.55<3.51): 전대역 기준 TBS(C≈24)의 TB-BLER 증폭이 실제 할당(~1/3 대역, C 적음)보다 비관적 → OLLA와 이중 보정. → 기준 할당을 "대표 per-UE 할당"으로 바꾸는 것이 개선 후보.


### 전이(transient) 진단 (1drop×1200slot, mbler+RBIR)

SE 3.51→**4.05**, BLER 0.224→**0.194** (500→1200슬롯). OLLA 전이가 일부 기여하나 수렴 후에도 BLER ~0.19로 타깃 0.1 대비 갭 잔존 → 남은 갭은 구조적: ① per-RB 입도(MATLAB은 per-subcarrier ESM), ② AWGN-가파른 커브가 per-RB SINR 모델 오차를 직격(레거시 테이블은 자체 생태계에 튜닝돼 있었음). 

**최종 권장**: `matlab_bler=0` 기본 유지. MATLAB 모드는 (a) OLLA step/window 재튜닝, (b) 임계값 기준할당 교정(대표 per-UE 할당), (c) subband-단위 ESM 입도 개선 후 재평가. RBIR 포팅 자체는 완료·검증됨(32/32)이며 남은 항목들의 전제 조건.


## 처리량 최대화 MCS 선택 (2026-07-04)

E.2 future-work 구현: `argmax over 28 MCS of TBS(m)×(1−TB_BLER(esinr,m,TBS(m)))`.

- **구현**: nr_l2sm `Build_TputMCS_Grid`(esinr −20~+40dB, 0.05dB 그리드 사전계산 → 스케줄러 조회 O(1)) + `determine_MCS` 두 복사본 훅. 플래그 `matlab_tput_mcs`(기본 0), 기준할당 `tput_mcs_ref_rbs`(0=auto num_rb/3 — (b)의 대표 per-UE 할당 반영).
- **그리드 검증**: 28개 MCS 전부 사용(CQI 맵은 15개), 경계 단조·~1dB 간격, 변조 전환 물리 타당(256QAM 시작 +19.8dB). **MCS27 경계 25.9dB ≈ 임계값 재생성의 독립 계산 25.86dB — 교차검증 일치.** 셀프테스트 448/448+32/32 유지.

### A/B (2drop×500slot)

| arm | SE | 1st-tx BLER |
|-----|----:|----:|
| legacy | 5.63 | 0.118 |
| RBIR + CQI맵 MCS | 3.26 | 0.228 |
| RBIR + 처리량최대 MCS | 3.34 | 0.201 |

**결과: 개선 미미(+2.7%, 노이즈 수준).** 원인 — **OLLA와 목적 충돌**: 처리량최대의 최적 동작점은 가파른 커브에서 BLER~0.2(무릎)인데 OLLA는 0.1을 향해 SINR 입력을 지속 하향 → 선택기가 최적점 아래 MCS로 밀림. 두 외루프가 상충. (OLLA-off 진단으로 순수 선택기 동작 확인 — 결과 아래 추가)


### OLLA-off 진단 — 충돌 가설 반박됨

tput+OLLA-off: SE **2.22** (OLLA-on 3.34 대비 **−34%**), BLER1tx 0.192.

**해석 교정**: OLLA가 최적점을 방해한다는 가설은 틀렸다. OLLA의 (음수) 오프셋은 오히려 **필수 보정**이었다 — 즉 진짜 지배 요인은 **스케줄러의 SINR 추정이 수신기의 실현 RBIR 유효 SINR보다 체계적으로 낙관적**이라는 것. OLLA를 끄면 낙관 추정 그대로 과공격 MCS 선택 → 재전송 소용돌이/드롭(1st-tx BLER에는 안 잡히는 손실) → SE 급락. 레거시 테이블은 커브가 완만해 이 추정 오차에 둔감했지만, MATLAB 커브(~1dB 워터폴)는 그대로 벌점을 매긴다.

**최종 결론 (MATLAB BLER 트랙 전체)**:
- 완성·검증된 것: 테이블 통합(448/448), RBIR ESM(32/32), 임계값 재생성, 처리량최대 MCS(그리드 교차검증) — 모두 스펙 재현 정확.
- 남은 지배 병목: **스케줄러 SINR 추정 ↔ 수신 실현 유효 SINR의 일관성**. 다음 레버는 MCS 격자/ESM이 아니라 추정 경로: (i) 수신기가 실측한 RBIR ESINR을 (지연) 피드백해 스케줄러 추정으로 사용, (ii) 또는 추정-실현 편향의 오프라인 캘리브레이션. OLLA의 단일 dB 오프셋은 이 보정의 조악한 근사로서 이미 그 방향으로 일하고 있음.
- **기본값 `matlab_bler=0` 유지.** 모든 플래그(rbir/tput_mcs/thresholds)는 연구용으로 보존.


## 실측 유효 SINR 피드백 (추정-일관성, 2026-07-04)

`matlab_esinr_fb`: 수신기의 실현 RBIR 유효 SINR과 스케줄러 추정(SCHEDULE_DECISION.capacity)의 차이를 UE별 적분 제어(gain 0.1, 클램프 ±[−15,+5]dB)로 누적해 링크어댑테이션 SINR 입력을 교정. OLLA 적용 3지점에 병기(역할 분담: 본 컨트롤러=추정편향 소거, OLLA=잔여 BLER 조절).

### 4-arm A/B (2drop×500slot)

| arm | SE | 1st-tx BLER |
|-----|----:|----:|
| legacy | 6.31 | 0.106 |
| matlab(RBIR) | 3.38 | 0.228 |
| + esinr_fb | 3.63 | 0.259 |
| **+ esinr_fb + tput_mcs** | **4.17** | 0.225 |

**판정**:
- **시너지 입증**: tput-MCS 단독은 +2.7%에 그쳤지만(입력 편향 탓), 편향 보정과 결합하니 **+23.5%** (3.38→4.17) — "촘촘한 MCS 래더는 입력이 정확할 때만 값을 한다"는 가설 확인. fb+tput > fb 단독(3.63) > tput 단독(3.34).
- MATLAB 모드 최적 스택 확정: `matlab_bler=1 + matlab_rbir=1 + matlab_esinr_fb=1 + matlab_tput_mcs=1`.
- **legacy 대비 갭 −34% 잔존**. 남은 후보: per-UE 스칼라 보정의 한계(RB/레이어별 편향 분포는 못 잡음), per-RB ESM 입도((c) 과제), 그리고 근본 질문 — legacy 테이블 자체가 관대(낙관)했을 가능성. 절대 진실 판정은 외부 참조(MATLAB SLS 예제 결과 등)와의 대조 필요.
- 기본값 전부 OFF 유지. 연구 스택으로 보존.


## OLLA 파라미터 최적화 시도 → 노이즈 벽 (2026-07-08)

MATLAB 스택(bler+rbir+esinr_fb+tput_mcs)에서 BLER이 목표 0.1 대비 ~0.22로 고착된 문제를 OLLA 재튜닝으로 풀려는 시도. **진단 계측 추가**: PerdropStatistics에 `initial_bler`, OLLA offset 평균/범위/양쪽 rail 포화 UE 수, esinr_fb 보정 분포 출력.

**구조적 발견 (config-무관, 견고)**: OLLA가 bang-bang으로 동작 — 60 UE 중 ~50개가 offset 상·하한 클램프에 포화, 중간 영역 ~10개뿐. 단일 per-UE dB 오프셋은 가파른 LDPC 워터폴 위에서 조악한 도구. **단, 실측 per-layer 1st-tx BLER은 모든 설정에서 ~0.10-0.13** — OLLA는 파라미터와 무관하게 평균 BLER 타깃을 이미 맞추고 있음(BLER 규제는 작동, SE 차이가 노이즈).

**스윕 (2drop×500slot, UMi 8port mx12 UE20)** — 라운드1(단독+5동시), 라운드2(4동시):
| arm | 변경 | SE | vs base |
|-----|------|----:|----:|
| base (w50,t0.1,mgn0.05,↓0.8/↑0.2,[−10,3]) | — | 4.076 | — |
| w100 (window100,mgn0.02) | 측정노이즈↓ | 4.386 | +7.6% |
| t15 (target0.15) | 목표완화 | 4.211 | +3.3% |
| gentle (↓0.2/↑0.05) | 스텝완화 | 4.151 | +1.8% |
| clamp20 (min−20) | 하한확장 | 4.045 | −0.8% |
| t20 (target0.2) | | 3.816 | −6.4% |

**노이즈 플로어 측정 (결정적)**: `R2_repro`는 base와 **완전 동일 cfg 재실행**인데 SE 4.076→3.858 (**−5.4%**), drop2만 4.135→3.674 (**−11%**). 원인: cfg가 100스레드 요청 × 다중 arm 동시 → 64코어 오버서브스크립션 → 알려진 RNG 레이스(PF 피드백이 스레드 인터리빙 비결정성 증폭, [[rng-race-correctness-not-reproducibility]]) 자극. **파라미터 효과(수%) < 노이즈(5.4%)** → 라운드1/2 순위 신뢰 불가.

**라운드3 결정판 (paired: 동일시드 4토폴로지, 20스레드×3 no-oversub, drop별 페어)**:
| drop | base | best(w100+t15) | combo(+mx8) |
|-----:|-----:|-----:|-----:|
| 1 | 4.285 | 4.482 | 4.539 |
| 2 | 3.588 | **1.843** | **5.087** |
| 3 | 3.537 | 4.406 | 3.459 |
| 4 | 3.513 | 3.799 | 3.083 |
| mean±SEM | 3.73±0.19 | 3.63±0.62 | 4.04±0.47 |
| BLER | 0.113 | 0.123 | 0.132 |

**drop2가 결정타**: 세 arm이 같은 시드로 **동일 토폴로지**를 봤는데 3.588/1.843/5.087로 2.7× 벌어짐. OLLA(완만한 링크어댑 보정)만으로 물리 불가 → **RNG 레이스가 시드 페어링마저 파괴**. 오버서브스크립션 제거(20스레드)로도 intra-run 스레드 레이스는 잔존(→ [[rng-race-correctness-not-reproducibility]]의 "thread-local fix로도 노이즈 불변"과 일치). SEM들이 전부 겹침 → **OLLA 파라미터 차이는 통계적으로 미해결**.

**최종 판정 (OLLA 트랙)**:
- **부정 결과(가치 있음)**: OLLA 파라미터 효과는 시뮬레이터 노이즈 플로어 아래. 4가지로 입증 — ①동일-cfg 재실행 5.4%, ②동일-토폴로지 2.7× 스윙, ③SEM 전부 중첩, ④실측 BLER은 파라미터 무관하게 타깃(~0.11). OLLA 튜닝은 SE 레버가 아님.
- **약한 best-estimate 권장** (통계 유의 아님, 점추정+물리적 근거): `olla_window_size 100 + olla_bler_margin 0.02` (윈도우↑ = 측정노이즈↓, 물리적으로 견고, 점추정 never-worse). 공격 방향(target 0.15 / max_offset 8)은 점추정 SE는 좋으나(combo mean 4.04 최고) BLER을 0.13으로 타깃 위로 밀고 분산이 큼 — 채택 시 주의. 현행 기본값 유지도 실질 동등.
- **진짜 레버는 OLLA가 아니라 재현성**: RNG 레이스의 per-work-item 시딩 수정(→ 결정론적 A/B 가능) 또는 drop 대량화(std~0.4-1.2 → 0.3 효과 해상에 ~60 drop 필요, 비현실적). 그 전까지 수% SE 비교는 원천적으로 불가.


## OLLA 최적화 결론 — 재현성 확보 후 반전 (2026-07-10)

RNG 레이스 3건(HARQ 코인·O2I 손실·간섭 빔 tilt)을 per-work-item 결정론 시딩으로 수정해
**완전 재현성**(동일 cfg 2회 bit-identical + 8스레드==64스레드 bit-identical)을 확보한 뒤,
노이즈 없는 클린 스윕(3drop×500slot, 동일시드 페어드, MATLAB 스택 bler+rbir+esinr_fb+tput_mcs)으로 재판정.

| arm | drop별 SE | 평균 | vs base | BLER(L1tx) | rail lo/hi |
|-----|-----------|-----:|--------:|-----------:|-----------:|
| base (w50,t0.1,[−10,+3]) | 4.037/3.593/4.118 | 3.916 | — | 0.103 | 22/28 |
| w100 (window100,mgn0.02) | 3.955/3.561/3.994 | 3.837 | −2.0% | 0.103 | 23/29 |
| t15 (target0.15) | 4.149/3.661/4.239 | 4.016 | +2.6% | 0.107 | 20/31 |
| combo (w100+t15+mx8) | 4.377/4.789/4.502 | 4.556 | +16.4% | 0.118 | 23/20 |
| **mx8 단독 (max_offset 8)** | **4.393/4.720/4.678** | **4.597** | **+17.4%** | 0.111 | 24/**18** |

**판정 (이전 "노이즈 벽" 결론의 반전)**:
- 노이즈 시대(±5.4% 재현 불가)에는 "OLLA 파라미터 효과 없음"이 결론이었으나, 결정론 확보 후
  **`olla_max_offset` 3→8 단독으로 +17.4% (3 drop 모두 승리, 페어드 일관)** — 진짜 레버였음.
- **메커니즘**: esinr_fb가 추정 편향을 소거하면 스케줄러 추정이 비관인 UE들이 다수 생기는데,
  이들의 OLLA 오프셋이 +3dB 상한에 포화(base 28/60 UE)되어 공격적 MCS 진입이 막혀 있었다.
  상한을 8dB로 열면 rail-hi 28→18로 풀리고 실측 BLER은 0.111로 목표 근처 유지.
- w100(−2.0%)/t15(+2.6%)는 mx8 위에 추가 이득 없음(combo 4.556 ≈ mx8 4.597) — 채택 불필요.

**최종 권장 (MATLAB 스택 전용)**: 기존 OLLA 설정에서 **`olla_max_offset 8.0` 한 줄만 변경**.
(window/target/step은 현행 유지. 레거시 BLER 모드에서는 미검증 — 적용 범위를 MATLAB 스택
cfg로 한정. 하한 −10은 비구속 확인됨.)

**부수 확정**: ① OLLA는 여전히 rail-분극 동작(중간영역 ~1/3)이나 평균 BLER 규제는 유효
② 스레드 수 불변성으로 A/B 프로토콜 단순화(스레드 수 무관) ③ 과거 legacy 6.31 vs MATLAB
4.17 갭 수치는 노이즈 시대 측정이라 재검증 필요 — mx8 반영 시 갭이 축소될 가능성
(legacy 재실행은 향후 과제).
