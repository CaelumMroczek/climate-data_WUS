import os
import s3fs
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray
import xesmf as xe
from tqdm import tqdm
import datetime

def extract_monthly_area_weighted_variable(
    variable,
    start_year,
    end_year,
    run,
    output_csv,
    coord_file,
    temp_dir,
    shapefile_path,
    buffer_m=500,
    exclude_basins=["11325000", "12143700"]
):
    os.makedirs(temp_dir, exist_ok=True)
    fs = s3fs.S3FileSystem(anon=True)

    if run == "hist":
        s3_prefix = "wrf-cmip6-noversioning/downscaled_products/gcm/cesm2_r11i1p1f1_historical/postprocess/d02/"
        file_template = f"{{var}}.daily.cesm2.r11i1p1f1.hist.d02.{{year}}.nc"
    else:
        s3_prefix = f"wrf-cmip6-noversioning/downscaled_products/gcm/cesm2_r11i1p1f1_ssp{run}/postprocess/d02/"
        file_template = f"{{var}}.daily.cesm2.r11i1p1f1.ssp{run}.d02.{{year}}.nc"

    coord_ds = xr.open_dataset(coord_file)
    lat2d = coord_ds["lat2d"].squeeze()
    lon2d = coord_ds["lon2d"].squeeze()

    lat_target = np.arange(np.floor(lat2d.min()), np.ceil(lat2d.max()), 0.05)
    lon_target = np.arange(np.floor(lon2d.min()), np.ceil(lon2d.max()), 0.05)
    lon_grid, lat_grid = np.meshgrid(lon_target, lat_target)

    target_grid = xr.Dataset({
        "lat": ("y", lat_grid[:, 0]),
        "lon": ("x", lon_grid[0, :])
    })

    gdf = gpd.read_file(shapefile_path)
    gdf = gdf[~gdf["basin_id"].isin(exclude_basins)]
    gdf = gdf.to_crs("EPSG:6933")
    gdf["geometry"] = gdf.geometry.buffer(buffer_m)
    gdf = gdf.to_crs("EPSG:4326")

    all_results = []

    for year in range(start_year, end_year + 1):
        if variable == "prec":
            files_to_download = [file_template.format(var=v, year=year) for v in ("prec", "t2max")]
        else:
            files_to_download = [file_template.format(var=variable, year=year)]

        local_files = {}

        for fname in files_to_download:
            s3_key = s3_prefix + fname
            local_path = os.path.join(temp_dir, fname)
            print(f"📅 Downloading {fname}")
            try:
                fs.get(s3_key, local_path)
                local_files[fname.split(".")[0]] = local_path
            except Exception as e:
                print(f"⚠️ Failed to download {fname}: {e}")

        try:
            if variable == "prec":
                prec_ds = xr.open_dataset(local_files["prec"])
                tmax_ds = xr.open_dataset(local_files["t2max"])

                prec = prec_ds["prec"]
                tmax = tmax_ds["t2max"] - 273.15  # Convert to Celsius

                time_coord = pd.to_datetime(prec_ds['day'].values.astype(int).astype(str), format="%Y%m%d")
                prec = prec.assign_coords(time=("day", time_coord))
                tmax = tmax.assign_coords(time=("day", time_coord))

                lat_da = xr.DataArray(lat2d.data, dims=("lat2d", "lon2d"), name="lat")
                lon_da = xr.DataArray(lon2d.data, dims=("lat2d", "lon2d"), name="lon")
                prec = prec.assign_coords(lat=lat_da, lon=lon_da)
                tmax = tmax.assign_coords(lat=lat_da, lon=lon_da)

                regridder = xe.Regridder(prec, target_grid, method="bilinear", periodic=False, reuse_weights=False)
                prec_ll = regridder(prec)
                tmax_ll = regridder(tmax)

                snow_daily = prec_ll.where(tmax_ll < 0.0)
                snow_monthly = snow_daily.groupby("time.month").sum()

                monthly = prec_ll.groupby("time.month").sum()

                for label, data in zip(["area_weighted_mean", "snow_prec"], [monthly, snow_monthly]):
                    data = data.rename({"x": "lon", "y": "lat"})
                    data = data.assign_coords({
                        "lon": ("lon", target_grid["lon"].values),
                        "lat": ("lat", target_grid["lat"].values)
                    })
                    data = data.rio.write_crs("EPSG:4326").rio.set_spatial_dims(x_dim="lon", y_dim="lat")
                    data_eq = data.rio.reproject("EPSG:6933")

                    res_m = abs(data_eq.rio.resolution()[0])
                    cell_area_km2 = (res_m ** 2) / 1e6

                    for _, row in gdf.iterrows():
                        basin_id = row.basin_id
                        geom_eq = gpd.GeoSeries([row.geometry], crs="EPSG:4326").to_crs("EPSG:6933").geometry[0]

                        for m in range(1, 13):
                            try:
                                clipped = data_eq.sel(month=m).rio.clip([geom_eq.__geo_interface__], data_eq.rio.crs, drop=True)
                                pixel_count = (~np.isnan(clipped)).sum().item()
                                if pixel_count == 0:
                                    raise ValueError("All NaN after clip")

                                total = np.nansum(clipped.values) * cell_area_km2
                                basin_area_km2 = pixel_count * cell_area_km2
                                area_weighted = total / basin_area_km2

                                result = {
                                    "basin_id": basin_id,
                                    "year": year,
                                    "month": m,
                                    "variable": variable,
                                    "stat": label,
                                    "value": area_weighted
                                }
                                all_results.append(result)
                            except Exception:
                                continue
            else:
                var_path = local_files[variable]
                data_ds = xr.open_dataset(var_path)
                var_data = data_ds[variable]

                if variable.startswith("t"):
                    var_data = var_data - 273.15  # Convert temperature to Celsius

                time_coord = pd.to_datetime(data_ds['day'].values.astype(int).astype(str), format="%Y%m%d")
                var_data = var_data.assign_coords(time=("day", time_coord))

                lat_da = xr.DataArray(lat2d.data, dims=("lat2d", "lon2d"), name="lat")
                lon_da = xr.DataArray(lon2d.data, dims=("lat2d", "lon2d"), name="lon")
                var_data = var_data.assign_coords(lat=lat_da, lon=lon_da)

                regridder = xe.Regridder(var_data, target_grid, method="bilinear", periodic=False, reuse_weights=False)
                var_ll = regridder(var_data)

                monthly = var_ll.groupby("time.month").mean()

                monthly = monthly.rename({"x": "lon", "y": "lat"})
                monthly = monthly.assign_coords({
                    "lon": ("lon", target_grid["lon"].values),
                    "lat": ("lat", target_grid["lat"].values)
                })
                monthly = monthly.rio.write_crs("EPSG:4326").rio.set_spatial_dims(x_dim="lon", y_dim="lat")
                monthly_eq = monthly.rio.reproject("EPSG:6933")

                res_m = abs(monthly_eq.rio.resolution()[0])
                cell_area_km2 = (res_m ** 2) / 1e6

                for _, row in gdf.iterrows():
                    basin_id = row.basin_id
                    geom_eq = gpd.GeoSeries([row.geometry], crs="EPSG:4326").to_crs("EPSG:6933").geometry[0]

                    for m in range(1, 13):
                        try:
                            clipped = monthly_eq.sel(month=m).rio.clip([geom_eq.__geo_interface__], monthly_eq.rio.crs, drop=True)
                            pixel_count = (~np.isnan(clipped)).sum().item()
                            if pixel_count == 0:
                                raise ValueError("All NaN after clip")

                            total = np.nansum(clipped.values) * cell_area_km2
                            basin_area_km2 = pixel_count * cell_area_km2
                            area_weighted = total / basin_area_km2

                            result = {
                                "basin_id": basin_id,
                                "year": year,
                                "month": m,
                                "variable": variable,
                                "stat": "area_weighted_mean",
                                "value": area_weighted
                            }
                            all_results.append(result)
                        except Exception:
                            continue

        except Exception as e:
            print(f"⚠️ Failed processing year {year}: {e}")
        finally:
            for f in local_files.values():
                os.remove(f)

    df = pd.DataFrame(all_results)

    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"✅ Done: {len(df)} rows saved to {output_csv}")

    return df
