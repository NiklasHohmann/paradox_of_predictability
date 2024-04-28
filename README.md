# paradox_of_predictability

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11059550.svg)](https://doi.org/10.5281/zenodo.11059550)

Supplementary code and data for "The paradox of Predictability Provides a Bridge Between 1 Micro- and Macroevolution".

## Authors (no specific order)

__Masahito Tsuboi__  
Lund University  
Web page: [www.masahitotsuboi.com](https://www.masahitotsuboi.com/)  
ORCID: [0000-0002-0144-2893](https://orcid.org/0000-0002-0144-2893)

__Niklas Hohmann__  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

__Kjetil L. Voje__  
ORCID:  [0000-0003-2556-3080](https://orcid.org/0000-0003-2556-3080)

__Melanie Hopkins__  
ORCID: [0000-0002-3580-2172](https://orcid.org/0000-0002-3580-2172)

__Meghan Balk__  
ORCID: [0000-0003-2699-3066](https://orcid.org/0000-0003-2699-3066)

__Sophie Nilén__  
ORCID: [0009-0002-0996-0182](https://orcid.org/0009-0002-0996-0182)

__Erik I. Svensson__  
ORCID: [0000-0001-9006-016X](https://orcid.org/0000-0001-9006-016X)

__Lee Hsiang Liow__  
ORCID: [0000-0002-3732-6069](https://orcid.org/0000-0002-3732-6069)

__Gene Hunt__  
ORCID: [0000-0001-6430-5020](https://orcid.org/0000-0001-6430-5020)

## Requirements

R (version >= 4) and the RStudio IDE

## Usage

Download the code from GitHub. In RStudio, open the file `paradox_of_predictability.Rproj`. This open the RProject, installs the `renv` package (if it is not already installed), and set the working directory correctly. Next, run

```r
renv::restore()
```

to install all required packages in the correct version. Now you are set up to interact with the code and data.

Then you can inspect and run the code in the directory `code/`. The main analysis is contained in the script `code/analysis.R`.

## Repository structure

* _code_ : folder with R code
* _data_ : folder with raw and modified data
  * contains various files with raw time series data, named after publication and year
  * _output_ : folder with data extracted from time series
  * _var_div_dat_ : folder with data extracted from time series, with added data on timescale (paleo or neo). This data was added manually
* _figs_ : folder for figures
* _renv_ : `renv` package folder
* _.gitignore_ : untracked files
* _.Rprofile_ : R session settings
* _LICENSE_ : Apache 2.0 license file
* _paradox_of_predicability.Rproj_ : R Project file
* _README.md_ : readme file
* _renv.lock_ : lock file for `renv` package

## License

Apache 2.0, see LICENSE file for details.

## Funding information

This work was funded by

1. The Center for Advanced Studies in Norwegian Academy of Science and Letters to Thomas F. Hansen and Christophe Pélabon and
2. The Japanese Society for the Promotion of Science, Overseas Research Fellowship to Masahito Tsuboi (Grant identification number: 202270014)
