import os
import s3fs
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray
import xesmf as xe
from tqdm import tqdm

def extract_area_weighted_variable(
    variable,
    start_year,
    end_year,
    run,
    output_csv,
    coord_file="/Users/caelum/Documents/GitHub/climate-data_downloads/wrfinput_d02_coord.nc",
    temp_dir="downloaded_nc_files",
    shapefile_path="/Users/caelum/Documents/GitHub/streamflow_WesternUS/1-data/shapefiles/headwater-catchments_shp/headwater_catchments.shp",
    buffer_m=500,
    exclude_basins=["11325000", "12143700"]
):
    os.makedirs(temp_dir, exist_ok=True)
    fs = s3fs.S3FileSystem(anon=True)

    # Determine S3 path prefix and filename format
    if run == "hist":
        s3_prefix = "wrf-cmip6-noversioning/downscaled_products/gcm/cesm2_r11i1p1f1_historical/postprocess/d02/"
        file_template = f"{variable}.daily.cesm2.r11i1p1f1.hist.d02.{{year}}.nc"
    else:
        s3_prefix = f"wrf-cmip6-noversioning/downscaled_products/gcm/cesm2_r11i1p1f1_ssp{run}/postprocess/d02/"
        file_template = f"{variable}.daily.cesm2.r11i1p1f1.ssp{run}.d02.{{year}}.nc"

    # Load WRF lat/lon
    coord_ds = xr.open_dataset(coord_file)
    lat2d = coord_ds["lat2d"].squeeze()
    lon2d = coord_ds["lon2d"].squeeze()

    lat_target = np.arange(np.floor(lat2d.min()), np.ceil(lat2d.max()), 0.05)
    lon_target = np.arange(np.floor(lon2d.min()), np.ceil(lon2d.max()), 0.05)
    lon_grid, lat_grid = np.meshgrid(lon_target, lat_target)

    target_grid = xr.Dataset({
        "lat": (["y", "x"], lat_grid),
        "lon": (["y", "x"], lon_grid)
    })

    # Load and buffer catchments
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf[~gdf["basin_id"].isin(exclude_basins)]
    gdf = gdf.to_crs("EPSG:6933")
    gdf["geometry"] = gdf.geometry.buffer(buffer_m)
    gdf = gdf.to_crs("EPSG:4326")

    all_results = []
    all_skipped = []

    for year in range(start_year, end_year + 1):
        filename = file_template.format(year=year)
        s3_key = s3_prefix + filename
        local_path = os.path.join(temp_dir, filename)

        print(f"📥 Downloading {filename}")
        try:
            fs.get(s3_key, local_path)
        except Exception as e:
            print(f"⚠️ Failed to download {filename}: {e}")
            continue

        try:
            data_ds = xr.open_dataset(local_path)
            var = data_ds[variable]

            # Create CF-compliant lat/lon 2D DataArrays
            lat_da = xr.DataArray(lat2d.data, dims=("lat2d", "lon2d"), name="lat")
            lon_da = xr.DataArray(lon2d.data, dims=("lat2d", "lon2d"), name="lon")
            var = var.assign_coords(lat=lat_da, lon=lon_da)

            # Regrid to regular lat/lon
            regridder = xe.Regridder(var, target_grid, method="bilinear", periodic=False, reuse_weights=False)
            var_ll = regridder(var)

            # Sum across time
            annual = var_ll.sum(dim="day", skipna=True)
            annual = annual.rename({"x": "lon", "y": "lat"})

            # Assign proper 1D lat/lon coordinates to match dimensions
            annual = annual.assign_coords({
                "lon": ("lon", target_grid["lon"].isel(y=0).values),
                "lat": ("lat", target_grid["lat"].isel(x=0).values)
            })

            # Add CRS and spatial dimensions
            annual = annual.rio.write_crs("EPSG:4326").rio.set_spatial_dims(x_dim="lon", y_dim="lat")
            annual_eq = annual.rio.reproject("EPSG:6933")

            res_m = abs(annual_eq.rio.resolution()[0])
            cell_area_km2 = (res_m ** 2) / 1e6

            for _, row in gdf.iterrows():
                shape_id = row.get("basin_id")
                geom = row.geometry

                try:
                    geom_eq = gpd.GeoSeries([geom], crs="EPSG:4326").to_crs("EPSG:6933").geometry[0]
                    clipped = annual_eq.rio.clip([geom_eq.__geo_interface__], annual_eq.rio.crs, drop=True)

                    if clipped.rio.width == 0 or clipped.rio.height == 0:
                        raise ValueError("No data in bounds.")

                    valid_mask = ~np.isnan(clipped.values)
                    pixel_count = valid_mask.sum()

                    if pixel_count == 0:
                        raise ValueError("All NaN after clip")

                    total_value = np.nansum(clipped.values) * cell_area_km2
                    basin_area_km2 = pixel_count * cell_area_km2
                    area_weighted_mean = total_value / basin_area_km2

                    all_results.append({
                        "basin_id": shape_id,
                        "year": year,
                        "variable": variable,
                        "area_weighted_mean": area_weighted_mean,
                        "valid_pixels": int(pixel_count),
                        "basin_area_km2": basin_area_km2
                    })

                except Exception as e:
                    all_skipped.append({"year": year, "basin_id": shape_id, "reason": str(e)})

        except Exception as e:
            print(f"⚠️ Failed processing {filename}: {e}")
        finally:
            os.remove(local_path)

    df = pd.DataFrame(all_results)
    df.to_csv(output_csv, index=False)
    print(f"✅ Done: {len(df)} entries saved to {output_csv} | Skipped: {len(all_skipped)}")
    return df
