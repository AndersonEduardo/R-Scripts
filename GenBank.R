##CODIGO PARA EXTRAIR SEQUENCIA DE NUCLEOTIDEOS NO GENBANK##

library("ape", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
cotton_acc <- c("U56806","U12712","U56810","U12732","U56802")
cotton <- read.GenBank(cotton_acc,species.names = T)
names_accs <- data.frame(species=attr(cotton,"species"),accs=names(cotton))
names(cotton) <- attr(cotton,"species")
write.dna(cotton,"cotton.fas",format="fasta") #nao funcionou...
