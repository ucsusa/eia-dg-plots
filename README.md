# Plot EIA data on distributed and centralized generation capacity
This code creates plots for an upcoming UCS blog series showing 
 - The comparitive growth of centralized and distributed generation
 - The growth of distributed generation by technology type
 - Comparison of installed cost of PV systems vs electricity prices for residential

## Sources

EIA, collected via the following PUDL data sets:
 - [Core EIA 861 Yearly distributed generation tech](https://data.catalyst.coop/pudl/core_eia861__yearly_distributed_generation_tech)
 - [Core EIA 861 Yearly net metering customer fuel class](https://data.catalyst.coop/pudl/core_eia861__yearly_net_metering_customer_fuel_class)
 - [Core EIA 861 Yearly non-net metering customer fuel class](https://data.catalyst.coop/pudl/core_eia861__yearly_non_net_metering_customer_fuel_class)
 - [Out EIA yearly generators](https://data.catalyst.coop/pudl/out_eia__yearly_generators)

EIA, collected via API:
 - [Average Price of Electricity to Ultimate Customers](https://www.eia.gov/opendata/browser/electricity/retail-sales?frequency=annual&data=price;&facets=sectorid;stateid;&sectorid=RES;&stateid=US;&sortColumn=period;&sortDirection=desc)

### Data that must be downloaded manually

LBNLs "Tracking the Sun" data is only available programmatically in a raw format in a datalake. To run this code, a subsest of the data must be downloaded manually:
1. Go to the [Tracking the Sun shiny app](https://emp-lbnl.shinyapps.io/tts_visualization/)
2. Go to the tab "Installed Prices"
3. Downbload the csv file for the first graph, "Installed Prices"
4. Save the csv to the directory `<path-to-repo>/workflow_solution/data/National_Residential_Installed Prices (2000-2023).csv`

## Using

The repo contains both a snamekmake workflow and Jupyter notebooks. The Jupyter notebooks produces additional plots.

## Dependencies

Both the notebook and snakemake workflow utilize a private repo to apply UCS styles to the final plots (`UCSmpl`). Two lines in the code will need to be removed to run the project without this library:

```python
import UCSmpl
...
with plt.style.context('ucs_light'):
```

After removing the second line be sure to adjust indentation accordingly.