## tutorial para trabalhar com Zonation + R
## FONTE: https://github.com/cbig/zonator
##        https://cbig.github.io/zonator/articles/zonator-project.html
##        https://cbig.github.io/zonator/articles/zonator-results.html

#install.packages("zonator")

# se instalar nao diretamente pelo install.packages nao funcionar, instalar pelo devtools:
#install.packages("devtools")
#library(devtools)
#devtools::install_github("cbig/zonator", ref = "develop", build_vignettes = TRUE)
#devtools::install_github("cbig/zdat")

library(zonator)
library(zdat)

# diretorio de trabalho
projectFolder = '/home/anderson/Projetos/Wil_zonation'

setwd(projectFolder)

# Define the names of the variants within the project
variant_names <- c("01_variant", "02_variant", "03_variant")

# Create a new project from scratch

#opcao 1: basico
create_zproject(name = "projeto_teste", dir = projectFolder, variants = variant_names)

#opcao 2: fornecendo diretorio com o raster da distribuicao das sps
input_raster_dir = '/home/anderson/Projetos/Wil - zonation/data'
new_project <- create_zproject(name = "projeto_teste", dir = projectFolder, 
                               variants = variant_names,
                               spp_template_dir = input_raster_dir,
                               spp_file_pattern = "^species[0-9].tif$")

#opcao 3: fornecendo pesos
input_raster_dir = '/home/anderson/Projetos/Wil_zonation/data'
new_project <- create_zproject(name = "projeto_teste", dir = projectFolder, 
                               variants = variant_names,
                               spp_template_dir = input_raster_dir,
                               spp_file_pattern = "^species[0-9].tif$",
                               weight = c(1, 1, 1, 2, 3, 2, 1))





# Start by creating a project using the tutorial data
#setup.dir <- system.file("extdata/basic", package = "zdat") #caminho pro deles
setup.dir <- paste(projectFolder, '/projeto_teste', sep='') #caminho pro meu
tutorial.project <- load_zproject(setup.dir)

# Get a specific variant, 01_core_area_zonation
variant.caz <- get_variant(tutorial.project, "01_variant")


# Let's also rename the features and groups while we're at it
featurenames(variant.caz) <- c("Koala", "Masked.owl", "Powerful.owl", 
                               "Tiger.quoll", "Sooty.owl", "Squirrel.glider",
                               "Yellow-bellied.glider")


