# paradox_of_predictability
Supplementary code and data for "The paradox of predictability"

## Authors

__Niklas Hohmann__  (creator and maintainer of repository)  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

__Masahito Tsuboi__  
Lund University  
Web page: [www.masahitotsuboi.com](https://www.masahitotsuboi.com/)  
ORCID: [0000-0002-0144-2893](https://orcid.org/0000-0002-0144-2893)

## Requirements

R (version >= 4) and the RStudio IDE

## Usage

Download the code from GitHub. In RStudio, open the file `paradox_of_predictability.Rproj`. This open the RProject, installs the `renv` package (if it is not already installed), and set the working directory correctly. Next, run

```r
renv::restore()
```

to install all required packages in the correct version. Now you are set up to interact with the code and data.

Then you can inspect and run the code in the directory `code/`

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

Add information here.
