Depth 발표 자료 — FLAG cut baseline d3c9272

발표 문서: 상위 폴더 presentation.txt
본문 7장: Sky depth 측정 → 20-band 2D map → 검출 천체 기반 depth와 비교.
m450 전체 이미지는 slide 4의 짧은 보조 화면, m475 막대그래프는 slide 7 질문 대응용이다.
현재 발표에 참조된 플롯만 보관한다.

Sky 방법 예시: COSMOS1 m650
- cosmos1_7dt_depth_m650_sky_mask_example.png: 동일 영역의 원본/마스킹 후 비교.
  baseline raw 중심과 크기 표를 사용하며, 원본 해상도의 셀 (행 4, 열 4)을 표시한다.
  prepare_sky_mask_example.py로 이 한 그림을 재현할 수 있다.
- cosmos1_7dt_depth_m650_aperture_flux_distribution.pdf: 기존 baseline flux 분포.

나머지 발표용 Depth 플롯
- cosmos1-4_20band_depth_maps_baseline.png
- source_depth_examples_baseline.png
- sky_source_offset_audit.png
- depth_comparison_by_band/m600_sky_source_depth_baseline.png
- depth_comparison_by_band/m575_sky_source_depth_baseline.png
- depth_comparison_by_band/m475_sky_source_depth_baseline.png

별도 COSMOS3_m450_ZP_check 폴더의 cosmos3_m450_full_image.jpg까지 총 9개 플롯이다.
plot_manifest.csv는 현재 보관 파일, 사용 슬라이드, 크기와 해시를 기록한다.
발표에 사용하지 않는 PNG/PDF 등은 삭제했다. 원본 output/ 아래의 플롯은 유지한다.
CSV, 분석 스크립트와 ipynb는 유지한다. 기존 분석 노트북을 전체 재실행하면
원래의 모든 밴드별 플롯이 다시 생성될 수 있으며, 현재 선택 목록과는 구분한다.
