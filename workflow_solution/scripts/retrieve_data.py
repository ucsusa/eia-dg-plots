import pandas as pd
from dotenv import load_dotenv
from gridstatus import EIA
import requests

load_dotenv("../.env")

if __name__ == "__main__":
        # Pull in the raw DG data from PUDL

        # Distributed generation is spread across three datasets, each showing
        # online capacity by year (i.e. not tracking additions and retirements)

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

        # Download the centralized generation data
        central_raw_url = (f"https://data.catalyst.coop/pudl/"
                f"out_eia__yearly_generators.csv?"
                f"_labels=on&_stream=on&generator_operating_date__notblank=1&"
                f"report_date__gte=2024-01-01&_size=max")

        central_raw = pd.read_csv(central_raw_url)
        central_raw.to_csv(snakemake.output.central_raw)

        # Get EIA data on average residential electric rates
        # https://www.eia.gov/opendata/browser/electricity/retail-sales?frequency=annual&data=price;&facets=sectorid;stateid;&sectorid=RES;&stateid=US;&sortColumn=period;&sortDirection=desc;

        # Instantiate the EIA client to get the API key
        eia = EIA()

        # API url
        url = (
                "https://api.eia.gov/v2/electricity/retail-sales/data/"
                "?frequency=annual"
                "&data[0]=price"
                "&facets[sectorid][]=RES"
                "&facets[stateid][]=US"
                "&sort[0][column]=period"
                "&sort[0][direction]=desc"
                "&offset=0"
                "&length=5000"
        )

        # Set up headers with the API key
        headers = {
                "X-Api-Key": eia.api_key
        }

        # Send the GET request to the API
        r = requests.get(url, headers=headers)
        r.raise_for_status()  # Raise an error if the request failed

        # Parse and load the response data into a DataFrame
        response = r.json()['response']
        print(f"Total records: {response['total']}")

        res_costs_raw = pd.DataFrame(response['data'])
        res_costs_raw.to_csv(snakemake.output.res_costs_raw)
