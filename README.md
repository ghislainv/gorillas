[![DOI](https://zenodo.org/badge/20937/ghislainv/gorillas.svg)](https://zenodo.org/badge/latestdoi/20937/ghislainv/gorillas)

# :gorilla: :congo_kinshasa: Estimating Grauer's gorilla distribution and population size in eastern Democratic Republic of the Congo

## Data

Data can be found in the `data` folder:

- The `rasters` folder includes a GeoTIFF raster stack file called `environ.tif` with all the environmental variables used for the species distribution model.
- The `gorillas` folder includes three .csv data files: 
    - `gorillas-erate-dens.csv`: data for the density/e-rate relationship
    - `gorillas-erate.csv`: data to compute local densities and the weighted mean density
    - `gorillas.csv`: presence-absence data for the species distribution model
    
## R script

The R script `gorillas.R` can be run to reproduce the results of the following study:

**Plumptre A. J., S. Nixon, D. Kujirakwinja, G. Vieilledent, R. Critchlow, E. A. Williamson, R. Nishuli, A. Kirkby and J. S. Hall**. 2016. Catastrophic Decline of World's Largest Primate: 80% Loss of Grauer's Gorilla (_Gorilla beringei graueri_) Population Justifies Critically Endangered Status. _PLoS One_. **11**(10): e0162697. [doi:[10.1371/journal.pone.0162697](https://doi.org/10.1371/journal.pone.0162697)].

To run the script, you can either execute the `gorillas.sh` shellscript file or run the R script under a R GUI such as RStudio.

## Results

All results are saved in the `results` and `report` folders. The `report.Rmd` file is used to embed R computations and text and report briefly the main results of the study in the file `report/report.pdf`. 

![Grauer's gorilla distribution](/results/SDA_ggmap.png)

## License

Data and R script are available under the GNU General Public License version 3 (see `LICENSE` file).
