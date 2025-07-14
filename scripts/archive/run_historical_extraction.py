from extract_precip_by_basin import extract_area_weighted_variable
import os
import pandas as pd

# Define historical range
start_year = 1980
end_year = 2014
years = range(start_year, end_year + 1)

# Define variables to extract
variables = ["prec", "t2", "t2min", "t2max"]

# Define run type
run_type = "hist"

# Output directory
output_dir = "outputs_historical_wide"
os.makedirs(output_dir, exist_ok=True)

# Loop through each variable
for var in variables:
    print(f"🔄 Processing variable: {var}")

    df_wide = pd.DataFrame()

    for year in years:
        print(f"   📆 Year: {year}")
        try:
            # Run the extraction for this year
            df = extract_area_weighted_variable(
                variable=var,
                start_year=year,
                end_year=year,
                run=run_type,
                output_csv=None  # disable writing per-year files
            )

            # Reshape to wide format: rows = years, columns = basin_ids
            df_pivot = df.pivot(index="year", columns="basin_id", values="area_weighted_mean")

            # Append this year's row to the wide dataframe
            df_wide = pd.concat([df_wide, df_pivot])

        except Exception as e:
            print(f"❌ Failed for {var} {year}: {e}")

    # Sort rows by year
    df_wide = df_wide.sort_index()

    # Save wide-format CSV
    output_csv = os.path.join(output_dir, f"{var}_hist_wide.csv")
    df_wide.to_csv(output_csv)
    print(f"✅ {var} saved to {output_csv}")
