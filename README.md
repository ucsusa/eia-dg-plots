# Plot EIA data on distributed and centralized generation capacity
This code creates plots for an upcoming UCS blog series showing 
 - The comparitive growth of centralized and distributed generation
 - The growth of distributed generation by technology type

## Sources

Data is from EIA, collected via the following PUDL data sets:
 - [Core EIA 861 Yearly distributed generation tech](https://data.catalyst.coop/pudl/core_eia861__yearly_distributed_generation_tech)
 - [Core EIA 861 Yearly net metering customer fuel class](https://data.catalyst.coop/pudl/core_eia861__yearly_net_metering_customer_fuel_class)
 - [Core EIA 861 Yearly non-net metering customer fuel class](https://data.catalyst.coop/pudl/core_eia861__yearly_non_net_metering_customer_fuel_class)
 - [Out EIA yearly generators](https://data.catalyst.coop/pudl/out_eia__yearly_generators)

## Using

The repo contains both a snamekmake workflow and a Jupyter notebook. The Jupyter notebook produces additional plots.

## Dependencies

Both the notebook and snakemake workflow utilize a private repo to apply UCS styles to the final plots (`UCSmpl`). Two lines in the code will need to be removed to run the project without this library:

```python
import UCSmpl
...
with plt.style.context('ucs_light'):
```
