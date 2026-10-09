# 7DT Single-band Matching 발표 계획

## 발표 목표

- 발표 길이: 약 5분
- 구성: 제목을 포함하여 총 5장
- 중심 질문: 서로 다른 band에서 검출된 source를 어떻게 동일 천체로 판정했고, 그 결과를 어떻게 검증했는가?
- 핵심 흐름: 20개 single-band catalog → 7DT master catalog → FINALCAT match → 검증된 union catalog

세부적인 conflict 처리, candidate fallback, catalog column 구성은 본 발표에서 제외한다.

---

## Slide 1. From Single-band Catalogs to a Master Catalog

### 목적

분석의 입력과 matching의 전체 구조를 한 번에 보여준다.

### 슬라이드 내용

- 20개 7DT band: `m400`–`m875`
- 각 band에서 독립적인 source detection 및 3-arcsec aperture photometry 수행
- quality selection과 valid-footprint selection 후 총 78,882개 band detection 사용
- matching은 두 단계로 수행

```text
20 single-band catalogs
        ↓
7DT cross-band matching
        ↓
19,489 7DT masters
        ↓
FINALCAT positional matching
        ↓
75,259 union masters
```

### 그림

다음 그림의 `Matching input count`와 `Adopted matching uncertainty` 패널만 crop하여 오른쪽에 작게 배치한다.

`../../output/COSMOS1/singlebandphot/matching/figures/2_singphotandmatching/singleband_input/singleband_input_summary.png`

### 발표 문장

> 각 band의 검출 수와 image quality가 다르기 때문에, 동일한 고정 반경이 아니라 band별 positional uncertainty를 반영하여 matching했습니다.

---

## Slide 2. 7DT Cross-band Matching Method

### 목적

7DT 내부 matching의 중심 원리만 설명한다.

### 슬라이드 내용

- 모든 band가 정렬된 common pixel grid에 있다고 가정
- `m400`부터 wavelength 순서로 master catalog를 점진적으로 구성
- 각 band의 positional uncertainty는 image FWHM과 S/N scale로 정의

$$
\sigma_b = \frac{\mathrm{FWHM}_b}{2.355\,\mathrm{S/N}_{\min,b}}
$$

- source와 현재 master의 combined uncertainty 사용

$$
\sigma_{\mathrm{match}} = \sqrt{\sigma_{\mathrm{master}}^2 + \sigma_b^2}
$$

- normalized separation이 3보다 작은 pair를 match candidate로 선택

$$
D = \frac{d}{\sigma_{\mathrm{match}}} < 3
$$

- matched source가 추가될 때 master position을 inverse-variance weighted mean으로 갱신

$$
x_{\mathrm{master}}
=
\frac{\sum_b x_b/\sigma_b^2}{\sum_b 1/\sigma_b^2}
$$

### 그림

PowerPoint에서 다음 정도의 간단한 schematic을 직접 제작한다.

```text
new-band source → 3σ candidate search → one-to-one match → weighted master update
```

### 발표 문장

> 핵심은 source separation을 band별 positional uncertainty로 정규화하고, 3-sigma 안에 있는 source만 같은 천체로 연결하는 것입니다.

---

## Slide 3. 7DT Matching Validation

### 목적

matching 결과가 최종 master position과 일관적인지 보여준다.

### 핵심 결과

| Metric | Result |
|---|---:|
| Input band detections | 78,882 |
| 7DT masters | 19,489 |
| Strict 20-band masters | 1,157 |
| Suspect members, normalized residual ≥ 3 | 117 |
| Suspect fraction | 0.148% |
| Median normalized residual | 0.364 |
| 95th percentile | 1.647 |
| 99th percentile | 2.333 |

각 accepted member와 최종 weighted master position 사이의 residual을 다음과 같이 검증한다.

$$
R = \frac{d_{\mathrm{member,final}}}
{\sqrt{\sigma_{\mathrm{final}}^2 + \sigma_b^2}}
$$

### 그림

다음 그림에서 두 패널만 크게 crop한다.

- `Assignment separation`
- `Final consistency`

`../../output/COSMOS1/singlebandphot/matching/figures/2_singphotandmatching/matching_7dt/7dt_matching_summary.png`

`Assignment separation`의 0 부근에는 초기 `m400` source가 포함되어 있으므로, 발표에서는 분포의 전체 형태와 cutoff에만 집중한다.

### 강조 문구

> 78,882개 accepted member 중 99.85%가 최종 master position과 3-sigma 이내에서 일치했습니다.

---

## Slide 4. FINALCAT Matching and Radius Validation

### 목적

FINALCAT match radius를 데이터 기반으로 선택했음을 보여준다.

### 슬라이드 내용

- 19,489개 7DT master와 valid footprint 내 58,747개 FINALCAT source를 sky coordinate에서 비교
- 주어진 radius 안의 pair를 separation 순으로 one-to-one assignment
- FINALCAT 위치를 10 또는 20 arcsec 이동한 8개 random-shift sample로 chance association 추정

| Radius | Actual matches | Shifted random | Random / actual |
|---:|---:|---:|---:|
| 0.5 arcsec | 2,682 | 55.8 ± 5.2 | 2.08% |
| **0.7 arcsec** | **2,977** | **107.9 ± 7.2** | **3.62%** |
| 1.0 arcsec | 3,212 | 211.4 ± 16.9 | 6.58% |

### 그림

아래 figure의 세 패널을 그대로 사용한다.

`../../output/COSMOS1/singlebandphot/matching/figures/2_singphotandmatching/finalcat_radius/finalcat_radius_diagnostics.png`

### 발표 문장

> 0.7 arcsec까지는 실제 match가 유의미하게 증가하지만, 1.0 arcsec에서는 추가 match보다 random association이 더 빠르게 증가합니다. 따라서 0.7 arcsec를 채택했습니다.

선택된 match의 separation은 median 0.205 arcsec, 95th percentile 0.578 arcsec이다.

---

## Slide 5. Final Results and Take-home Message

### 최종 결과

| Master class | Number |
|---|---:|
| 7DT + FINALCAT matched | 2,977 |
| 7DT only | 16,512 |
| FINALCAT only | 55,770 |
| Total union catalog | 75,259 |

### 그림

다음 그림에서 `Sky distribution`과 `Master classes` 패널만 사용한다.

`../../output/COSMOS1/singlebandphot/matching/figures/2_singphotandmatching/master_population/master_population_summary.png`

### Take-home messages

1. 7DT 내부 matching은 band별 positional uncertainty를 반영한 normalized 3-sigma criterion을 사용했다.
2. 내부 consistency 검증에서 accepted member의 99.85%가 3-sigma 이내였다.
3. Random-shift test를 통해 FINALCAT match radius를 0.7 arcsec로 선택했다.
4. 최종적으로 2,977개의 7DT–FINALCAT counterpart와 75,259개의 union master catalog를 얻었다.

### 마지막 발표 문장

> Positional uncertainty와 random-association test를 함께 사용하여, 내부적으로 일관되고 외부 catalog contamination을 통제한 master catalog를 구축했습니다.

---

## 제작 원칙

- 한 슬라이드에 핵심 숫자는 최대 3개만 색으로 강조한다.
- 기존 multi-panel figure를 축소해서 통째로 넣지 말고 필요한 panel만 crop한다.
- 방법론에서는 식을 모두 읽지 않고 `band-dependent uncertainty`, `normalized 3-sigma`, `weighted master` 세 표현만 강조한다.
- 결과에서는 catalog 규모보다 `99.85% internal consistency`와 `0.7 arcsec radius validation`을 가장 크게 보여준다.
- candidate conflict, magnitude-based ambiguity ranking, fallback status 등은 질문이 있을 때만 설명한다.
