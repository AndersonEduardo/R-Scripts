##fonte: http://rcompanion.org/handbook/G_11.html

install.packages('ordinal')
install.packages('lsmeans')

library(ordinal)
library(lsmeans)

Input =("
Instructor  Question       Likert
 Fuu         Informative     8
 Fuu         Informative     9
 Fuu         Informative     9
 Fuu         Informative     8
 Fuu         Delivery       10
 Fuu         Delivery        9
 Fuu         Delivery        8
 Fuu         Delivery        8
 Fuu         VisualAides     7
 Fuu         VisualAides     7
 Fuu         VisualAides     6
 Fuu         VisualAides     7
 Fuu         AnswerQuest     8
 Fuu         AnswerQuest     9
 Fuu         AnswerQuest     9
 Fuu         AnswerQuest     8
 Jin         Informative     7
 Jin         Informative     8
 Jin         Informative     6
 Jin         Informative     5
 Jin         Delivery        8
 Jin         Delivery        8
 Jin         Delivery        9
 Jin         Delivery        6
 Jin         VisualAides     5
 Jin         VisualAides     6
 Jin         VisualAides     7
 Jin         VisualAides     7
 Jin         AnswerQuest     8
 Jin         AnswerQuest     7
 Jin         AnswerQuest     6
 Jin         AnswerQuest     6
 Mugen       Informative     5
 Mugen       Informative     4
 Mugen       Informative     3
 Mugen       Informative     4
 Mugen       Delivery        8
 Mugen       Delivery        9
 Mugen       Delivery        8
 Mugen       Delivery        7
 Mugen       VisualAides     5
 Mugen       VisualAides     4
 Mugen       VisualAides     4
 Mugen       VisualAides     5
 Mugen       AnswerQuest     6
 Mugen       AnswerQuest     7
 Mugen       AnswerQuest     6
 Mugen       AnswerQuest     7
")

Data = read.table(textConnection(Input),header=TRUE)

### Order levels of the factor; otherwise R will alphabetize them
Data$Instructor = factor(Data$Instructor,
                         levels=unique(Data$Instructor))

### Create a new variable which is the likert scores as an ordered factor
Data$Likert.f = factor(Data$Likert,
                       ordered = TRUE)

### modelo estatistico
model = clm(Likert.f ~ Instructor + Question + Instructor:Question,
            data = Data,
            threshold="symmetric")

### teste de significancia dos parametros do modelo
anova(model,
      type = "II")

### teste posteriori (para comparar entre os grupos)
marginal = lsmeans(model,
                      pairwise ~ Instructor + Question,
                      adjust="tukey")        ### Tukey-adjusted comparisons

marginal$contrasts




