if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import pandas as pd
    import UCSmpl

    df_combined = pd.read_csv(snakemake.input.combined_data)
    df_distributed_tech = pd.read_csv(snakemake.input.distributed_tech_data)

    # DG vs centralized plot
    start_year = 1900
    df_plot = df_combined[(df_combined['year'] >= start_year) & (df_combined['year'] <= 2023)]

    with plt.style.context('ucs_light'):
        plt.figure(figsize=(10, 6))    
        plt.plot(df_plot['year'], df_plot['central_mw'], label='Centralized capacity (MW)')
        plt.plot(df_plot['year'], df_plot['dg_capacity_mw'], label='Distributed capacity (MW)')
        plt.xlabel('Year')
        plt.ylabel('Capacity [MW]')
        plt.title(f"Centralized and Distributed Capacity ({start_year}–2023)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('log')
        plt.savefig(snakemake.output.dg_vs_central)

    # DG plot by type
    df_plot = df_distributed_tech[df_distributed_tech['year'] >= 2010].copy()

    # Combine tech classes and clean up the names for plotting
    rename_map = {
        'all_storage': 'Storage',
        'chp_cogen': 'Fossil Fuel',
        'combustion_turbine': 'Fossil Fuel',
        'fuel_cell': 'Fossil Fuel',
        'hydro': 'Other',
        'internal_combustion': 'Fossil Fuel',
        'pv': 'Solar',
        'steam': 'Fossil Fuel',
        'wind': 'Other',
        'storage_nonpv': 'Storage',
        'storage_pv': 'Storage',
        'virtual_pv': 'Solar',
        'other' : 'Other'
    }
    df_plot['tech_class'] = df_plot['tech_class'].replace(rename_map)
    df_plot = df_plot.dropna(subset=['tech_class', 'dg_capacity_mw'])
    df_plot = df_plot.pivot_table(index='year', columns='tech_class', values='dg_capacity_mw', aggfunc='sum')
    df_plot = df_plot.drop(columns=['total'], errors='ignore')
    desired_order = ['Solar', 'Fossil Fuel', 'Storage', 'Other']
    df_plot = df_plot[desired_order]

    # Plot
    with plt.style.context('ucs_light'):
        ax = df_plot.plot(figsize=(10, 6))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))
        plt.ylabel('Capacity (MW)')
        plt.xlabel('Year')
        plt.title('Distributed Generation Capacity by Technology (2010-2023)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(snakemake.output.dg_by_type)
