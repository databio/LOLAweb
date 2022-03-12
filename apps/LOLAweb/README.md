## LOLAweb

### Getting Started

`LOLAweb` is an interactive application that implements functions from the LOLA package in a web interface via the [Shiny](https://shiny.rstudio.com/) framework.

The app is hosted and available for use at:

<http://lolaweb.databio.org>

To run the application locally, you'll need to first clone this repository:

```
git clone https://github.com/databio/LOLAweb.git
```

The app requires R to be installed, as well as several packages. From within R run the following to install the dependencies:

```
install.packages(c("ggplot2", "shiny", "DT", "shinyWidgets", "shinyjs", "sodium", "devtools", "shinyBS"))
devtools::install_github("databio/simpleCache")
devtools::install_github("databio/GenomicDistributions")

source("https://bioconductor.org/biocLite.R")
biocLite(c("LOLA", "GenomicRanges"))
```

You'll also need underlying data that is not available in this repository in order to establish the universes, example user sets and reference genome directories. Detailed guidance for how to download and organize reference data is available in the [LOLAweb project README](https://github.com/databio/lolaweb/#readme).

Note that one additional step for running the app locally willbe to either 1) pass the `$LWREF` and `$LWLOCAL` environment variables to your R session or 2) overwite the definition of `localDir` and `refDir` in `misc.R` with the respective paths on your machine.

With dependencies installed and data configuration complete you can now launch the app with `shiny::runApp()` from within R at the root of this directory.