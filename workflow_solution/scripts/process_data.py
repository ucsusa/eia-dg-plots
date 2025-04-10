import pandas as pd

# Function for data cleaning PUDL. Sums capacity by year using the date in the
# indicated column.
def annual_totals(df, date_field = 'report_date'):

    # Convert 'report_date' to datetime and extract year
    df[date_field] = pd.to_datetime(df[date_field])
    df['year'] = df[date_field].dt.year

    # Delete empty rows
    df = df.loc[df['capacity_mw'] != 0]

    # Group by 'report_year' and calculate total capacity
    result = df.groupby('year')['capacity_mw'].sum().reset_index()

    return result

if __name__ == "__main__":

    frames = []
    for file in [snakemake.input.dg_raw, snakemake.input.nm_raw, snakemake.input.nnm_raw]:
        df = pd.read_csv(file)
        frames.append(annual_totals(df))

    df_distributed = pd.concat(frames)

    # Group by 'year' and sum the 'capacity' values. This df now has annual online capacity
    df_distributed = df_distributed.groupby('year', as_index=False)['capacity_mw'].sum()

    # Rename capacity
    df_distributed = df_distributed.rename(columns={'capacity_mw': 'dg_capacity_mw'})

    central_raw = pd.read_csv(snakemake.input.central_raw)
    additions = annual_totals(central_raw, "generator_operating_date")
    retirements = annual_totals(central_raw, "generator_retirement_date")
    additions = additions.rename(columns={'capacity_mw': 'additions_mw'})
    retirements = retirements.rename(columns={'capacity_mw': 'retirements_mw'})
    df_central = pd.merge(additions, retirements, on='year', how='outer').fillna(0)
    df_central['net_change_mw'] = df_central['additions_mw'] - df_central['retirements_mw']
    df_central['central_mw'] = df_central['net_change_mw'].cumsum()
    
    # Combine the two data sets to plot
    df_combined = pd.merge(df_distributed, df_central, on='year', how='outer').fillna(0)

    df_combined.to_csv(snakemake.output.processed_data)