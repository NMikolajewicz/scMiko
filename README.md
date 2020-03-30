### scMiko

Set of general and module-specific utility functions used for scRNAseq analysis. Builds upon Seurat Package. 

## Installation

Firstly, please install or update the package devtools by running

```
install.packages("devtools")
```

Then the scMiko can be installed via

```
library(devtools)
devtools::install_github(
    repo = "NMikolajewicz/scMiko",
    ref = "master",
    auth_token = "a3c1c9b15c496991c952d1fe3ccc52db770f22fa")

private token: a3c1c9b15c496991c952d1fe3ccc52db770f22fa
```
## Built With

* [Seurat](https://satijalab.org/seurat/) - Seurat R toolkit for single cell genomics
* [flexdashboard](https://rmarkdown.rstudio.com/flexdashboard/) - Interactive dashboards for R

## Authors

* [Nicholas Mikolajewicz](https://scholar.google.ca/citations?user=LBWQMXsAAAAJ&hl=en&oi=ao)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hong Han and Kevin Brown from the [Moffat Lab](http://moffatlab.ccbr.utoronto.ca/)
