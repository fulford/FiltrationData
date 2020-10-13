# FiltrationData

This MATLAB script generates synthetic data for the Filtration 2 experiment in the 3rd year chemical engineering lab.
This data is exported to a Microsoft Excel file.

If you happen to be a 3rd year chemical engineering student at the Unversity of Edinburgh, you should note that your error analysis should **not** cover details revealed below (directly).

## Background

This script assumes that most data successfully follows the trends expected from the equations.
Small deviations are created as expected from experimental data:

- Random fluctuations are included in each time reported.
- A number of the inital points may be calculated for pressures other than the set pressure.
- Each run uses a concentration that features a deviation from the input calcium carbonate slurry concentration.
- A unique medium and cake resistance are used.

## Warnings

No back-calculation is performed on the resulting dataset - you should self sanity check that the values obtained are realistic.

The inital presssures are not restricted to being increasing (or decreasing) this can lead to some interesting inital points on the curves - to allow evaluation of this the curves *will* be automatically plotted folloing the data export.

The Excel sheet in this repostiory must be present as the location for MATLAB to save the values to - the cell locations are hard-coded.
