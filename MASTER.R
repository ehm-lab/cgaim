library(devtools)
library(usethis)
library(testthat)

package_name <- "cgaim"
package_path <- sprintf("C:/Users/PierreMasselot/Documents/Recherche/# R/Packages/%s",
  package_name)

# Package creation (run once)
# create_package(package_path)

# Set the working directory in the package folder
# setwd(package_path)

# Load all functions from package
load_all()

# Create tests
use_testthat()
use_test("basic_models", open = F)
use_test("build_constraint")
use_test("bootstrap")
use_test("confint")

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
usethis::use_package("mgcv")
usethis::use_package("gratia")
usethis::use_package("methods")
usethis::use_package("doParallel")
usethis::use_package("foreach")

# CMD check
devtools::check()

# Quick installation
devtools::install(quick = T)

# Package installation
remove.packages("cgaim")
devtools::install_github("PierreMasselot/cgaim")
library(cgaim)
