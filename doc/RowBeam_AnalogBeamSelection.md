# Row/Column-Beam Allocation — 아날로그 빔 선정 알고리즘

협력기관 제안 알고리즘의 구현 문서. subarray-partitioned 하이브리드 어레이에서
**매 drop마다 UE들의 RSRP 분포 투표로 수직 포트 행(row)별 아날로그 빔을 배정**한다.
cfg 스위치 `row_beam_enable`로 on/off (off = 기존 동작 bit-identical).

## 1. 배경 및 동기

기존 시뮬레이터는 아날로그 빔을 **per-링크**로 모델링한다: 각 UE가 자신에게 최적인
`(sector_z, sector_a)` 빔을 선택하고, 그 UE의 유효 채널(W_tx)은 전 포트에 그 빔을
적용한다. 그러나 실제 subarray-partitioned 패널은 **한 시점에 패널당 한 벌의 아날로그
가중치**만 걸 수 있으므로, MU-MIMO로 여러 UE를 같은 slot에 서빙할 때 UE마다 다른
아날로그 빔을 주는 것은 물리적으로 불가능한 낙관 가정이다.

이 기능은 성능 "향상" 레버가 아니라 **모델 교정**이다:

- 기대 방향: (ii) 제안 방식 < (i) 기존 per-링크 모델 (현실화로 SE 하락)
- 성공 기준: (ii) 제안 방식 > (iii) 전 행 단일 빔 (행 다양성의 실질 기여 확인)

## 2. 알고리즘 (drop마다, 섹터별)

수직 포트 행 V = `BS_Mp`개, 수평 포트 열 H = `BS_Np`개. 포트 (mi, ni)는
서브어레이 소자 (mi·K+k, ni·L+l), K = `BS_M/BS_Mp`, L = `BS_N/BS_Np`를 묶는다.

### ① UE 투표
섹터에 접속한 각 UE는 자기 기준 슬라이스(자신의 best azimuth `sector_a`, 자신의
UE빔 `(z,a,p)`)에서 zenith 빔별 RSRP를 보고, **best ± X dB 이내** 빔들에
RSRP 내림차순으로 **최대 K개** 투표한다 (tie는 빔 인덱스 오름차순).

- X = `row_beam_x_db` (기본 3.0 dB), K = `row_beam_max_cand` (기본 3)
- RSRP는 이미 빔스캔(`find_best_tx_beam`)에서 전량 캐시된
  `signal_RSRP_gain[sec][z][a][ue_z][ue_a][p]` 룩업 — 새 물리 계산 없음

### ② 득표 집계 → ③ 비례 배분
섹터 득표 벡터 count_z[]를 largest-remainder 방식으로 V개 행에 배분한다:
floor(득표비율×V) 우선, 잔여 슬롯은 소수부 내림차순(tie: 득표수 내림차순 → 인덱스
오름차순). 접속 UE가 없는 섹터는 boresight(zenith 90°) 최근접 빔에 전량 폴백.

### ④ 배치
배분된 빔을 **zenith 내림차순(가장 아래로 향하는 빔부터) 최상단 행부터** 연속 배치.
소자 z좌표가 mi에 비례(d_tx.z = m·dV)하므로 최상단 행 = mi = V−1. 같은 빔이 인접
행에 붙어 eType II 수직 기저와의 부정합을 최소화한다. 결과는
`sector[flat].row_beam_z[mi]`에 저장.

### 수평(azimuth) 처리 — `row_beam_az_mode`
- **0** (기본): 열은 기존 per-UE azimuth 유지 — 수직 배정 효과만 격리 측정
- **1**: 전 열 boresight 고정 (그리드에서 0 rad 최근접, tie는 낮은 인덱스)
- **2**: azimuth에도 ①~③ 동일 투표/배분 → 열 `sector[flat].col_beam_a[ni]`에
  azimuth 오름차순 배치 (2D 일반형; 4×4/4×8 등 서브어레이 실험 대비)

### 일관성 갱신
배정 후 각 UE의 서빙 RSRP를 **혼합 평균**으로 재계산한다. 포트별 RSRP 게인은 그
포트 빔만의 함수이므로 어레이 전체 값은 정확히 산술평균이다:

```
gain_mix = (1/V) Σ_mi gain[z_row[mi]][a*]          (az_mode 0/1)
gain_mix = (1/(V·H)) Σ_mi Σ_ni gain[z_row[mi]][a_col[ni]]   (az_mode 2)
```

`link_RSRP` / `link_antgain` / `str_signal` / `static_gain`을 갱신해 이후의
`Get_CouplingLoss`(coupling loss, geometry, SIR)와 RI 계산이 자동으로 일관된다.
셀 접속(select_serving_cell)은 best-빔 기준을 유지한다 (브로드캐스트 빔 접속 근사).

### 간섭 반영 (사용자 결정 ②)
- **광대역 간섭**: DU/Rural 경로에서 간섭 셀의 랜덤 zenith 빔 룩업을 그 셀의 **실제
  배정 빔 행 평균**으로 교체 (`Get_CouplingLoss`). 랜덤 draw 자체는 유지되어 RNG
  스트림이 flag on/off 무관하게 동일 — bit-identity 보장의 핵심.
- **시간영역 채널**: `Fourier_Transform_WithBF`의 W_tx가 간섭 링크(adj_sector ≠
  서빙)에도 그 섹터의 배정 빔을 사용하므로 간섭 채널 계수에 자동 반영.

### Zenith 그리드 오버라이드 (사용자 결정 ③)
`row_beam_num_zenith` > 0이면 `[min, max]°` 균등 B개 빔으로 교체 (azimuth 그리드는
불변). 예: 60–120° B=4 → 60/80/100/120°.

시나리오별 기본 그리드 주의:
- **Dense_Urban_eMBB_A (Config A)**: full-sim 그리드가 **1×1 고정빔**(az 0°, zen 102°)
  — 오버라이드 없이는 row-beam이 no-op (전 행 동일 빔 = 기존과 동일, 검증에 활용).
- **DU Config B/C, Rural**: 8az × 4zen (zenith 22.5/67.5/112.5/157.5°) — UMi에서
  22.5°/157.5°는 사실상 낭비(실효 2빔)라 오버라이드 권장.

## 3. cfg 파라미터

| 파라미터 | 기본값 | 설명 |
|---|---|---|
| `row_beam_enable` | 0 | 마스터 스위치. 0 = 기존 동작 bit-identical |
| `row_beam_az_mode` | 0 | 0=per-UE azimuth 유지, 1=boresight 고정, 2=열별 투표 배분 |
| `row_beam_force_uniform` | 0 | 1 = 전 행/열에 population argmax 빔 (A/B arm iii) |
| `row_beam_max_cand` | 3 | UE당 최대 투표 빔 수 K |
| `row_beam_x_db` | 3.0 | 투표 윈도우 X dB |
| `row_beam_zenith_min_deg` | 0 | 그리드 오버라이드 하한 (deg) |
| `row_beam_zenith_max_deg` | 0 | 그리드 오버라이드 상한 (deg) |
| `row_beam_num_zenith` | 0 | 오버라이드 빔 수 B (0 = 기존 그리드 유지) |

가드: `Calibration_mode==1` 또는 InH(TYPE==11) 조합은 시작 시 에러 종료.
주의: cfg 주석에 파라미터명 토큰을 쓰지 말 것 (`Get_parameter`의 토큰 스캔 특성).

## 4. 구현 위치

| 파일 | 내용 |
|---|---|
| `src/Initiallization.cpp` | `Assign_Row_Beams()` 신규 (투표/배분/배치/일관성 갱신, `Link_configuration` 끝에서 직렬 호출); zenith 그리드 오버라이드 + boresight 인덱스 (`Generate_bs_2D_DFT_beam_precoder`); `Get_CouplingLoss` DU 분기 간섭 룩업 |
| `src/channel_update.cpp` | `Fourier_Transform_WithBF` W_tx 빌드 — 행별 `row_beam_z[mi]`, (az_mode 2) 열별 `col_beam_a[ni]` |
| `src/channel.cpp` | dead `Update_v2`의 동등 블록 2곳 일관성 유지 (현재 호출자 없음) |
| `src/setSimulParam.cpp` | cfg 파싱, 가드, 시작 echo |
| `h/sector.h` | `row_beam_z[16]` / `col_beam_a[16]` 멤버 |
| `h/const.h`, `h/channel.h`, `h/Initiallization.h`, `h/common.h` | zenith 차원 4→8 확장 (`MAX_RSRP_SEC_ZENITH`, `signal_RSRP_gain`, `w/v/virtualization_weight_wv`, tilt 배열) + extern 크기 불일치 버그 수정 (`bs_tilt_*_LCS`) |
| `src/PerdropStatistics.cpp` | 기존 버그 수정: Schedule_map 진단 합산이 미스케줄 슬롯(`ue_selected == NO_UE(-1)`)으로 `links[-1]`을 OOB 읽던 문제 — 센티널 스킵. 이 수정으로 out.txt의 Schedule_map 4개 진단 값만 기존 대비 달라짐 (물리 출력은 전부 bit-identical) |

결정론: `Assign_Row_Beams`는 난수를 쓰지 않고 UE를 인덱스 순으로 스캔한다
(OMP 루프에서 무동기화로 채워지는 `sector[].ue_in_control`은 순서 비결정이라 사용 금지).

## 5. 검증

1. **flag-off bit-identity**: 수정 전/후 바이너리, 동일 cfg·시드 → 전체 출력 diff 동일
2. **flag-on sanity**: Σalloc=V 체크(위반 시 즉시 종료), 2회 실행 diff 동일(결정론),
   단일 UE + `force_uniform` + `max_cand 1`에서 배정 = 그 UE의 best 빔 전 행이고
   혼합 RSRP = 기존값 → coupling loss CDF가 flag-off와 동일 (end-to-end 등식)
3. **3-arm A/B** (동일 시드): (i) enable=0 / (ii) enable=1 az_mode=0 /
   (iii) enable=1 force_uniform=1. 지표: 평균/5%ile UE SE, 셀 SE, Wideband SIR·
   geometry CDF, 섹터별 배정 히스토그램. 판정: (ii)−(iii)가 헤드라인(행 다양성 기여),
   (i)−(ii)는 낙관 모델 교정폭.

## 6. A/B 결과

### 스모크 검증 (구현 시점, 짧은 런: 1 drop × 12 slot, 7 site, 2 UE/sector, DU-A + zenith 60–120° B=4)

| 검증 항목 | 결과 |
|---|---|
| flag-off bit-identity (수정 전 vs 후 바이너리) | ✔ 모든 물리 출력(CDF/throughput/geometry/SIR/CP loss) 동일. out.txt의 Schedule_map 진단 4개 값만 변경 — `links[-1]` OOB 수정의 의도된 효과 |
| flag-on 결정론 (2회 실행) | ✔ 전 출력 동일 |
| 퇴화 no-op 등식 (Config A 1×1 그리드, flag-on) | ✔ 전 행 동일 빔 → 물리 출력이 flag-off와 완전 동일 (혼합 RSRP·간섭 평균 경로 end-to-end 검증) |
| Σalloc = V 체크 | ✔ 위반 없음 |
| 2-drop 재배정 | ✔ drop마다 투표/배정 갱신 확인 |
| force_uniform / az_mode 2 | ✔ 단일 빔 전 행 / cols 배분 동작 확인 |
| 방향성 | 스케줄 UE SINR: (i) 14.47 dB > (ii) 투표 14.12 dB (낙관 교정) · 커버리지 median geometry: (ii) 35 dB > (iii) 단일빔 32 dB, SIR 25 > 22 dB (행 다양성 기여) |

### 본 A/B (장기 런)

(실행 후 추가 — 3-arm: (i) enable=0 / (ii) az_mode=0 / (iii) force_uniform=1, 동일 시드)

## 7. 알려진 한계 / 리스크

- eType II 코드북의 수직 기저(N2 DFT)는 행 간 균일 응답을 가정 — 행별 상이 빔이
  이를 깨서 양자화 손실이 늘 수 있다. 이는 arm (ii)가 측정하는 효과의 일부이며,
  연속 행 배치가 최선의 완화다.
- 접속은 best-빔 RSRP 기준 유지 — 배정 빔 기준 재접속(핸드오버 모델)은 범위 외.
- InH(TYPE==11) 및 Calibration 모드는 지원하지 않음 (시작 가드).
