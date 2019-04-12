library(googleComputeEngineR)

## setup
project = 'instancia-r'
zone = 'southamerica-east1-a'
account_key = '/home/anderson/Projetos/SDM megafauna Sul-Americana/instancia-R-7c898a713ed7.json'
Sys.setenv(GCE_AUTH_FILE = account_key,
           GCE_DEFAULT_PROJECT_ID = project,
           GCE_DEFAULT_ZONE = zone)
options(googleAuthR.scopes.selected = "https://www.googleapis.com/auth/cloud-platform")
gce_auth()


## global project
gce_global_project(project)
gce_global_zone(zone)
gce_get_project('instancia-r')

## instanciando a maquina virtual
vm = gce_vm(template = 'rstudio',
            name = 'rstudio-vm',
            username = 'AndersonEduardo',
            password = 'senha123123',
            predefined_type = 'n1-highmem-2')

gce_vm('rstudio-vm') #verificando

## finalizando job da mv
job = gce_vm_stop('rstudio-vm')

## listando  e inicializando instancias criadas previamente
gce_list_instances()
gce_vm_start('rstudio-vm')