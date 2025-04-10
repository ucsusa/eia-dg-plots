import pandas as pd

if __name__ == "__main__":
    ## Pull in the raw DG data from PUDL

    # Distributed generation is spread across three datasets, each showing total
    # online capacity per year (i.e. does not track additions and retirements)

    # Pre-2016 DG
    dg_raw_url = (f"https://data.catalyst.coop/pudl/"
            f"core_eia861__yearly_distributed_generation_tech.csv?"
            f"_labels=on&_stream=on&_size=max")

    dg_raw = pd.read_csv(dg_raw_url)
    dg_raw.to_csv(snakemake.output.dg_raw)

    # Net-metering data
    nm_raw_url = (f"https://data.catalyst.coop/pudl/"
            f"core_eia861__yearly_net_metering_customer_fuel_class.csv?"
            f"_labels=on&_stream=on&_size=max")

    nm_raw = pd.read_csv(nm_raw_url)
    nm_raw.to_csv(snakemake.output.nm_raw)


    # Non-net metering data
    nnm_raw_url = (f"https://data.catalyst.coop/pudl/"
            f"core_eia861__yearly_non_net_metering_customer_fuel_class.csv?"
            f"_labels=on&_stream=on&_size=max")

    nnm_raw = pd.read_csv(nnm_raw_url)
    nnm_raw.to_csv(snakemake.output.nnm_raw)


    ## Download the centralized generation data
    central_raw_url = (f"https://data.catalyst.coop/pudl/"
            f"out_eia__yearly_generators.csv?"
            f"_labels=on&_stream=on&generator_operating_date__notblank=1&"
            f"report_date__gte=2024-01-01&_size=max")

    central_raw = pd.read_csv(central_raw_url)
    central_raw.to_csv(snakemake.output.central_raw)
