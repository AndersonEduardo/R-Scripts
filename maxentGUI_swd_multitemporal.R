## Breve tutorial para MaxEnt usando dados multitemporais


## MODO 1: MaxEnt GUI (interface com botoes) usando SWD (Sample With Data) ##


# Fase 1: preparando dados de ocorrencia no R

#abrindo planilha de dados, ja limpos e organizados. OBS: colunas devem ser necessariamente: especie, longitude, latitude, idade do fossil
dataSetRaw = read.csv(file='/home/anderson/Downloads/teste_lipa/dataset_teste2.csv', header=TRUE, dec='.', sep=',')

dataSet = dataSetRaw[, c('Species', 'Longitude', 'Latitude', 'age')] #pegando as colunas certas, quando for o caso
names(dataSet) = c('species', 'lon', 'lat', 'age') #padronizando os nomes das colunas, de acordo com o necessario/interesse

#verificacao
str(dataSet)

# referencia:
# a coluna 'species' deve ser do tipo Factor; 
# a coluna 'lon' deve ser do tipo num
# a coluna 'lat' deve ser do tipo num
# a coluna 'age' deve ser do tipo num ou tipo int
# se nao, precisa abrir a planilha (preferencialmente com o
# bloco de notas) e procurar os ajustes necessarios.

#obtando dados das variaveis ambientais nos pontos de ocorrencia
source('/home/anderson/R-Scripts/paleoextract.R') #carregando funcao necessaria

dataSet$id = 1:nrow(dataSet) #por enquanto, a funcao paleoestract precisa de um 'id' para cada ponto para poder funcionar  : )
dataSetFull = paleoextract(x = dataSet, cols = c('lon', 'lat', 'age'), path = '/home/anderson/Downloads/teste_lipa/variaveis')
dataSetFull$species = 'Bradypus variegatus' #retornado a coluna que foi apagada pelo paleoextract
dataSetFull$id = NULL #apagando a coluna acessoria 'id'

#organizando a ordem das colunas (IMPORTANTE: precisa ser a mesma ordem das colunas do arquivo de background points, a ser criado na sequencia)
dataSetFull = dataSetFull[, c('species', 'lon', 'lat', 'age', sort(grep(pattern='bio', x=names(dataSetFull),value=TRUE)) )]

#salvando no computador
write.csv(x = dataSetFull, file = '/home/anderson/Downloads/teste_lipa/preguica_gigante/ocorrencias.csv', row.names = FALSE)



# Fase 2: backgroud points

#carregando funcao necessaria
source('/home/anderson/R-Scripts/paleobg.R')

#gerando background points
bgDataSet = paleobg(x = dataSet, colNames = names(dataSet), envFolder = '/home/anderson/Downloads/teste_lipa/variaveis', n = 10000)

#gerando o nome (necessario para o MaxnEnt GUI)
bgDataSet$species = 'background'

#organizando a ordem das colunas (IMPORTANTE: precisa ser a mesma ordem das colunas do arquivo de ocorrencias, criado ateriormente)
bgDataSetFull = dataSetFull[, c('species', 'lon', 'lat', 'age', sort(grep(pattern='bio', x=names(dataSetFull),value=TRUE)) )]

#salvando no computador
write.csv(x = bgDataSetFull, file = '/home/anderson/Downloads/teste_lipa/preguica_gigante/background.csv', row.names = FALSE)


# Fase 3: rodando no MaxEnt GUI

# agora, abra o MaxEnt GUI
# na janela para ' samples file': carregar o arquivo 'ocorrencias.csv'
# na jenela 'environmental layers': carregar o arquivo 'background.csv'
# ajustar a configuracao do MaxEnt de acordo com o interesse, normamente
# rodar e tchau!! :)


