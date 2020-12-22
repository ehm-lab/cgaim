library(devtools)
library(usethis)
library(testthat)

package_name <- "cgaim"
# package_path <- sprintf("C:/Users/masselpl/Documents/Recherche/# R packages/%s", 
#   package_name)

# Package creation (run once)
# create_package(package_path)

# Set the working directory in the package folder
# setwd(package_path)

# Load all functions from package
load_all()

# Create tests
use_testthat()
use_test("basic_models", open = F)

test()

# Vignettes
# use_vignette("SmoothMethod_simulation")

# Compile documentation
document()

# Adding dependencies
usethis::use_package("graphics")
usethis::use_package("stats")
usethis::use_package("scam")
usethis::use_package("scar")
usethis::use_package("cgam")
usethis::use_package("quadprog")
usethis::use_package("osqp")
usethis::use_package("limSolve")
usethis::use_package("Matrix")
usethis::use_package("grDevices")
usethis::use_package("graphics")
usethis::use_package("methods")
usethis::use_package("stats")
usethis::use_package("MASS")
usethis::use_package("cgam")

# CMD check
devtools::check()

# Package installation
install_github("PierreMasselot/cgaim")
library(cgaim)
