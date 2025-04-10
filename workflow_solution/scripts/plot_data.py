if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import pandas as pd
    import UCSmpl

    df_combined = pd.read_csv(snakemake.input.processed_data)
 
    # Filter for years
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
        plt.savefig('log-scale.png', dpi=100)
        plt.show()

    plt.savefig(snakemake.output.plot)