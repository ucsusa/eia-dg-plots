if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    from matplotlib.ticker import MaxNLocator
    import pandas as pd
    import UCSmpl

    df_combined = pd.read_csv(snakemake.input.combined_data)
    df_distributed_tech = pd.read_csv(snakemake.input.distributed_tech_data)
    df_costs = pd.read_csv(snakemake.input.costs_data)

    # DG vs centralized plot
    start_year = 1900
    df_plot = df_combined[(df_combined['year'] >= start_year) &
                          (df_combined['year'] <= 2023)]

    with plt.style.context('ucs_light'):
        plt.figure(figsize=(10, 6))    
        plt.plot(df_plot['year'], df_plot['central_mw'],
                 label='Centralized capacity (MW)')
        plt.plot(df_plot['year'], df_plot['dg_capacity_mw'],
                 label='Distributed capacity (MW)')
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
        'other': 'Other'
    }
    df_plot['tech_class'] = df_plot['tech_class'].replace(rename_map)
    df_plot = df_plot.dropna(subset=['tech_class', 'dg_capacity_mw'])
    df_plot = df_plot.pivot_table(index='year', columns='tech_class',
                                  values='dg_capacity_mw', aggfunc='sum')
    df_plot = df_plot.drop(columns=['total'], errors='ignore')
    desired_order = ['Solar', 'Fossil Fuel', 'Storage', 'Other']
    df_plot = df_plot[desired_order]

    # Plot
    with plt.style.context('ucs_light'):
        ax = df_plot.plot(figsize=(10, 6))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x,
                                                           _: f'{int(x):,}'))
        plt.ylabel('Capacity (MW)')
        plt.xlabel('Year')
        plt.title('Distributed Generation Capacity by Technology (2010-2023)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(snakemake.output.dg_by_type)

    # Costs plot
    with plt.style.context('ucs_light'):
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax2 = ax1.twinx()

        # Define colors
        pv_color = "#3044B5"
        price_color = "#CB2C30"

        # PV cost on the left axis
        ax1.plot(df_costs["Year"], df_costs["PV_Median_Cost"], color=pv_color,
                 label="Installed PV Cost")
        ax1.set_xlabel("Year")
        ax1.set_ylabel("Installed PV Cost ($/W)", color=pv_color)
        ax1.tick_params(axis="y", labelcolor=pv_color)
        ax1.set_ylim(0, 16)
        ax1.grid(True)

        # Set x-axis limits based on Year, and ensure it's only integer values
        ax1.set_xlim(df_costs["Year"].min()-1, df_costs["Year"].max()+1)
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

        # Electric cost on the right axis  
        ax2.plot(df_costs["Year"], df_costs["Residential_Electricity_Price"],
                 color=price_color, label="Electricity Price")
        ax2.set_ylabel("Electricity Price (cents/kWh)", color=price_color)
        ax2.tick_params(axis="y", labelcolor=price_color)
        ax2.set_ylim(0, 16)

        plt.title("Residential Installed PV Cost " +
                  "vs Electricity Price (2001-2023)")
        fig.tight_layout()
        plt.savefig(snakemake.output.pv_vs_elec_cost)
