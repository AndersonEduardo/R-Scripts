## tutorial para trabalhar com Zonation + R
## FONTE: https://github.com/cbig/zonator
##        https://cbig.github.io/zonator/articles/zonator-project.html
##        https://cbig.github.io/zonator/articles/zonator-results.html

#install.packages("devtools")
#library(devtools)
#devtools::install_github("cbig/zonator", ref = "develop", build_vignettes = TRUE) #instala o zonator
#devtools::install_github("cbig/zdat") #instala o zdat

library(zonator)
library(zdat)

# diretorio de trabalho
projectFolder = '/home/anderson/Projetos/Wil_zonation'

setwd(projectFolder)

# Define the names of the variants within the project
variant_names <- c("01_variant", "02_variant", "03_variant") #tipo "replicas" da analise - precisamos ver melhor este item

# Create a new project from scratch

#opcao 1: basico
create_zproject(name = "projeto_teste", dir = projectFolder, variants = variant_names)

#opcao 2: fornecendo diretorio com o raster da distribuicao das sps
input_raster_dir = '/home/anderson/Projetos/Wil - zonation/data' #ou seja, diretorio com os arquivos da distribuicao das especies
new_project <- create_zproject(name = "projeto_teste", dir = projectFolder, 
                               variants = variant_names,
                               spp_template_dir = input_raster_dir,
                               spp_file_pattern = "^species[0-9].tif$")

#opcao 3: fornecendo pesos
input_raster_dir = '/home/anderson/Projetos/Wil_zonation/data' #ou seja, diretorio com os arquivos da distribuicao das especies
new_project <- create_zproject(name = "projeto_teste", dir = projectFolder, 
                               variants = variant_names,
                               spp_template_dir = input_raster_dir,
                               spp_file_pattern = "^species[0-9].tif$",
                               weight = c(1, 1, 1, 2, 3, 2, 1))

## FEITO TUDO ISSO: é só abrir o Zonation e rodar o projeto gerado! : )
