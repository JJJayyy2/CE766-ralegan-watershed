# CE 766 Watershed Management: Evolutionary Optimization of Intervention Portfolios
## Ralegan Siddhi Watershed, Parner, Ahmednagar, Maharashtra

**Course:** CE 766 Watershed Management (Autumn 2026)  
**Instructor:** Riddhi Singh  
**Project Type:** Preliminary Results Submission (5 marks)  
**Submission Date:** October 22, 2026

---

## Overview

Multi-objective optimization study of watershed intervention portfolios in Ralegan Siddhi (Ahmednagar, Maharashtra). Uses NSGA-II to identify cost-effective combinations of check dams, farm ponds, contour bunds, and afforestation that maximize groundwater recharge while minimizing soil erosion over 2000–2023.

**Key finding:** 120 Pareto-optimal portfolios span ₹1–24 M cost, 382–501 TCM/yr recharge (+3% to +35% vs baseline), and 3,120–14,900 t/yr sediment (−81% to +10%).

---

## Files

### Data
- `rain.csv` — Daily gridded rainfall (IMD, 2000–2023, 8,766 rows)
- `pareto_portfolios.csv` — 120 Pareto-optimal solutions with decision variables and objectives

### Code
- `ralegan_portfolio_opt.py` — NSGA-II optimization script (Python 3.8+)

### Figures
- `preliminary_results.png` — 3-panel Pareto trade-off plot (300 dpi)

---

## Quick Start

```bash
pip install pymoo numpy pandas matplotlib xarray netcdf4
python ralegan_portfolio_opt.py rain.csv
```

Outputs: `pareto_portfolios.csv`, `preliminary_results.png`

---

## Results Summary

| Baseline | Pareto (mean) | Pareto (best recharge) |
|----------|---------------|------------------------|
| Recharge: 371.6 TCM/yr | 437.9 TCM/yr (+18%) | 500.9 TCM/yr (+35%) |
| Sediment: 16,394 t/yr | 8,880 t/yr (−46%) | 3,120 t/yr (−81%) |
| Cost: ₹0 | ₹12.1 M | ₹23.8 M |

### Intervention Selection (% of Pareto solutions)

- **Contour bunds:** 78–100% (universal)
- **Farm ponds:** 70–98% (mid/lower zones)
- **Check dams:** 60–72% (upper zones)
- **Afforestation:** 22–47% (lower cost-effectiveness)

---

## Methodology

**Hydrological model:** SCS-CN runoff + water-balance recharge (check dams, farm ponds, bunds) + RUSLE erosion

**Optimization:** NSGA-II, population 120, 150 generations, budget cap ₹25 M

**Objectives:** Minimize cost & sediment; Maximize recharge & irrigated area

---

## Data Sources

- IMD gridded rainfall: https://imdweb.lwcc.in/
- Future calibration: WRIS (https://indiawris.gov.in/wris/), CGWB, Bhuvan, ISRIC

---

## Limitations

- Parameters from literature (needs local calibration)
- No well-level validation yet (planned)
- Simplified livelihood model
- Assumes no structural degradation

---

## Next Steps

1. Acquire well-level data → calibrate recharge (WTF method)
2. High-res DEM + LULC → refine CN, K, LS, C
3. Scenario discovery & sensitivity analysis
4. Compare optimized vs. actual portfolios

---

## References

Deb et al. (2002) NSGA-II; Renard et al. (1997) RUSLE; USDA (1972) CN method; Gupta et al. (2009) KGE; Beven & Binley (1992) GLUE
