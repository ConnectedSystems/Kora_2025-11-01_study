# Paper dataset summary

Numbers derived from the fully-filtered pipeline CSVs
(`offshore_north/overall/` and `torres_strait/overall/`).
Filters applied on top of `ecorrap_expanded.parquet`:
1. `lin_ext > 0` — excludes partial-mortality / shrinking colonies (growth only)
2. Taxa filter — only taxa mapping to one of the 5 TARGET_GROUPS functional groups retained

---

## Offshore North

| Statistic | Value |
|-----------|-------|
| Reefs | Lizard, Moore |
| Transitions (survey years) | 2021_2022, 2022_2023 |
| Growth observations | 6772 |
| Survival observations | 10652 |
| Monitoring sites (reef × habitat_type) | 7 |
| Plots (reef × habitat_type × plot) | 43 |
| Taxa — growth | 33 |
| Taxa — survival | 34 |
| Habitat types | 4 (back, flank, front, lagoon) |
| Mean survey interval — growth (days) | 320.7 (SD: 86.3; range: 118–488) |
| Mean survey interval — survival (days) | 329.4 (SD: 75.3; range: 118–488) |


---

## Torres Strait

| Statistic | Value |
|-----------|-------|
| Reefs | Aukane, Dungeness, Masig |
| Transitions (survey years) | 2021_2022, 2022_2023 |
| Growth observations | 3524 |
| Survival observations | 5839 |
| Monitoring sites (reef × habitat_type) | 10 |
| Plots (reef × habitat_type × plot) | 51 |
| Taxa — growth | 30 |
| Taxa — survival | 31 |
| Habitat types | 4 (back, flank, front, lagoon) |
| Mean survey interval — growth (days) | 336.5 (SD: 5.6; range: 323–357) |
| Mean survey interval — survival (days) | 337.0 (SD: 5.7; range: 323–357) |


---

## Observation counts by reef and transition


### Offshore North — growth

| reef | transition | n |
|------|------------|---|
| lizard | 2021_2022 | 1543 |
| lizard | 2022_2023 | 1711 |
| moore | 2021_2022 | 1830 |
| moore | 2022_2023 | 1688 |
| **total** | | **6772** |

### Offshore North — survival

| reef | transition | n |
|------|------------|---|
| lizard | 2021_2022 | 1850 |
| lizard | 2022_2023 | 3515 |
| moore | 2021_2022 | 2219 |
| moore | 2022_2023 | 3068 |
| **total** | | **10652** |

### Torres Strait — growth

| reef | transition | n |
|------|------------|---|
| aukane | 2021_2022 | 674 |
| aukane | 2022_2023 | 680 |
| dungeness | 2021_2022 | 365 |
| dungeness | 2022_2023 | 113 |
| masig | 2021_2022 | 947 |
| masig | 2022_2023 | 745 |
| **total** | | **3524** |

### Torres Strait — survival

| reef | transition | n |
|------|------------|---|
| aukane | 2021_2022 | 875 |
| aukane | 2022_2023 | 1365 |
| dungeness | 2021_2022 | 480 |
| dungeness | 2022_2023 | 165 |
| masig | 2021_2022 | 1203 |
| masig | 2022_2023 | 1751 |
| **total** | | **5839** |

---

## Study variable summary (pipeline-filtered)

Sourced from the pipeline CSVs (post lin_ext and taxa filters).
Diameter ranges from `diam` (start-of-interval size, for both growth and survival obs.).
Depth: both regions use plot-level survey depth (`depth_cont`) from the IPM photogrammetry
file, now available for all reefs including Masig. Torres Strait also has continuous logger
mean depth (`depth_min_mean`) from EcoRRAP in-situ loggers (available at Masig Reef only) for
comparison.
Temperature: mean of daily maxima for ON (`temp_max_mean`), mean of daily means for TS
(`temp_mean_mean`). Colony ID is only tracked in the Torres Strait photogrammetry dataset.

| Identifier | Description | Offshore North<br>Count / [Range] | Torres Strait<br>Count / [Range] |
|------------|-------------|-----------------------------------|----------------------------------|
| `diam` | Estimated coral diameter at start of observation interval (cm), for both growth and survival obs. | Growth: 6772 / [0.5 – 131.83]<br>Survival: 10652 / [0.2 – 131.83] | Growth: 3524 / [0.2 – 190.45]<br>Survival: 5839 / [0.2 – 190.45] |
| `plot` | Physical monitoring plot (reef × site × habitat × depth × plot/quadrat label) | 43 | 51 |
| `depth` | Continuous depth (m): ON uses plot-level survey depth (`depth_cont`) from IPM photogrammetry file; TS uses logger mean depth (`depth_min_mean`) from EcoRRAP in-situ loggers (Masig Reef only) | 16 / [2.9 – 10.27] | 2 / [11.93 – 12.03] |
| `temp` | Mean of daily temperature (°C) from EcoRRAP in-situ loggers; daily-max mean for ON, daily mean for TS | 13 / [NaN – NaN] | 11 / [28.38 – 28.71] |
| `habitat_type` | Reef exposure category (`habitat_type`, e.g. "back"/"front") | 4 (back; flank; front; lagoon) | 4 (back; flank; front; lagoon) |
| `lon` | Longitude (decimal degrees) | 16 / [145.44 – 146.25] | 22 / [142.91 – 143.46] |
| `lat` | Latitude (decimal degrees) | 16 / [-16.88 – -14.65] | 21 / [-10.04 – -9.73] |
| `site` | Reef-level site code | 2 | 3 |
| `reef` | Monitored reef | 2 (Lizard; Moore) | 3 (Aukane; Dungeness; Masig) |
| `taxa` | Individual taxonomic group | Growth: 33<br>Survival: 34 | Growth: 30<br>Survival: 31 |
| `functional_group` | Morphological grouping of taxa | 5 | 5 |
| `colony_id` | Individual coral identifier | N/A | Growth: 2527<br>Survival: 3945 |

