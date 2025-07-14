
# WUS-D3 Climate Data Extraction Scripts

These scripts are designed to extract **monthly, area-weighted climate variables** (including precipitation, snow precipitation, and temperatures) from the WUS-D3 WRF-downscaled CMIP6 CESM2 datasets. The data is aggregated by hydrologic basin and suitable for ecohydrologic trend analyses, hydrologic modeling, or climate impact studies.

---

##  Script 1: `extract_monthly_variables.py`

### Purpose:
This script defines the function `extract_monthly_area_weighted_variable()` to:
- Download daily climate NetCDF files from the WUS-D3 S3 bucket (AWS).
- Regrid curvilinear WRF data to a regular lat/lon grid.
- Clip and aggregate values to buffered HUC8 basin polygons.
- Output monthly values of area-weighted means by basin.

### Key Features:
- **Supports multiple climate variables**: `prec`, `t2max`, `t2min`, `tmean`, etc.
- **Special handling for precipitation**:
  - Calculates **snow precipitation** as daily precipitation where daily max temperature < 0°C.
  - Outputs both total precipitation and snow precipitation.
- **Temperature values are converted to Celsius** automatically.
- Output includes `basin_id`, `year`, `month`, `variable`, `stat`, and `value`.

### Required Inputs:
- `variable`: Name of the climate variable (e.g., `"prec"`, `"t2max"`).
- `start_year`, `end_year`: Year range to extract.
- `run`: `"hist"` or a scenario run name (`"ssp245"`, `"ssp585"`, etc.).
- `output_csv`: Output file path (or `None` to return only DataFrame).
- `coord_file`: Local path to a NetCDF file containing WRF lat2d/lon2d coordinates.
- `temp_dir`: Temporary local directory for downloaded files.
- `shapefile_path`: Path to a shapefile containing polygons with a `basin_id` field.

---

##  Script 2: `run_monthly-climate_historical.py` (example batch runner)

### Purpose:
This script loops through multiple years and/or variables, calling `extract_monthly_area_weighted_variable()` for each combination and saving outputs to CSV.

### Example Logic:
```python
variables = ["prec", "t2max", "t2min"]
start_year = 1980
end_year = 1982
run_type = "hist"

for var in variables:
    df_all = pd.DataFrame()
    for year in range(start_year, end_year + 1):
        df = extract_monthly_area_weighted_variable(
            variable=var,
            start_year=year,
            end_year=year,
            run=run_type,
            ...
        )
        df_all = pd.concat([df_all, df], ignore_index=True)

    output_csv = os.path.join(output_dir, f"{var}_monthly_area_weighted.csv")
    df_all.to_csv(output_csv, index=False)
```

---

## Output Format

Each output CSV contains the following columns:

| basin_id | year | month | variable | stat               | value  |
|----------|------|--------|----------|--------------------|--------|
| 15020001 | 1980 | 1      | prec     | area_weighted_mean | 23.12  |
| 15020001 | 1980 | 1      | prec     | snow_prec          | 4.88   |
| 15020001 | 1980 | 1      | t2max    | area_weighted_mean | 4.17   |

---

## Notes & Considerations

- The script assumes **WRF coordinates** are provided separately via a file (not embedded in each climate file).
- For `prec`, both `prec` and `t2max` are downloaded to compute snow.
- Temperature units are automatically converted from **Kelvin to Celsius**.
- Basins are buffered by default (`500 m`) before clipping to avoid edge artifacts.
- Output is projected using **EPSG:6933** (Equal Area) for accurate area-weighting.

---

##  Dependencies
Make sure the following packages are installed in your environment:
```bash
xarray, rioxarray, geopandas, pandas, numpy, s3fs, xesmf, tqdm
```

---
