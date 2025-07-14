#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 13 11:42:47 2025

@author: caelum
"""

# initial install of packages
pip install s3fs xarray numpy pandas geopandas rioxarray xesmf tqdm

import os
import pandas as pd
from extract_monthly_variables import extract_monthly_area_weighted_variable

# Configuration
start_year = 1980
end_year = 2013
variables = ["prec", "t2", "t2min", "t2max"]
run_type = "hist"
coord_file = "/Users/caelum/Documents/GitHub/climate-data_WUS/wrfinput_d02_coord.nc"
temp_dir = "downloaded_nc_files"
shapefile_path = "/Users/caelum/Documents/GitHub/streamflow_WesternUS/1-data/shapefiles/headwater-catchments_shp/headwater_catchments.shp"
output_dir = "outputs_monthly_test"
os.makedirs(output_dir, exist_ok=True)

# Run batch extraction
for var in variables:
    print(f' Processing variable: {var}')
    df_all = pd.DataFrame()

    for year in range(start_year, end_year + 1):
        print(f'    Year: {year}')
        try:
            df = extract_monthly_area_weighted_variable(
                variable=var,
                start_year=year,
                end_year=year,
                run=run_type,
                output_csv=None,
                coord_file=coord_file,
                temp_dir=temp_dir,
                shapefile_path=shapefile_path
            )
            df_all = pd.concat([df_all, df], ignore_index=True)
        except Exception as e:
            print(f"❌ Error processing {var} {year}: {e}")

    # Save result
    out_csv = os.path.join(output_dir, f"{var}_monthly_hist.csv")
    df_all.to_csv(out_csv, index=False)
    print(f'✅ Saved {var} to {out_csv}')
