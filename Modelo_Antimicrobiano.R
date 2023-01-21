# remove.packages(c("ggplot2", "data.table"))
# install.packages('ggplot2', dep = TRUE)
# install.packages('data.table', dep = TRUE)
# install.packages("reticulate")
# install.packages("caret")
# install.packages("dplyr")
# install.packages("Tplyr")
# install.packages("readr")
# install.packages("tidyverse")
# install.packages("caret")


## Importando librerias

library(pROC)
library(lattice)
library(caret)
library(gmodels)
library(caret)
library(Tplyr)
library(plyr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(doMC)
library(recipes)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
# library(reticulate)
library(kernlab)
library(LiblineaR)
# library(Formula)
# library(truncnorm)
# library(brnn)
# library(monmlp)

registerDoMC(cores = 8)

## Para grabar el output en un fichero

output = 'output_antimicrobiano.txt'
try(output <- file(output, open="wt")) 
sink(output, type = "output", append = TRUE, split = TRUE)
sink(output, type = "message", append = TRUE)
source("output_antimicrobiano.txt", echo=TRUE, max.deparse.length=10000)

## Importando los archivos csv de descriptores
Descriptores2D <- read.table("./DescriptorsCSV/Descriptores2D_Dataframe.csv", 
                             sep = ",", 
                             fill = TRUE, 
                             header = TRUE)
Descriptores3D <- read.table("./DescriptorsCSV/Descriptores3D_Dataframe.csv", 
                             sep = ",", 
                             fill = TRUE, 
                             header = TRUE)
DescriptoresLipinski <- read.table("./DescriptorsCSV/DescriptoresLipinski_Dataframe.csv", 
                                   sep = ",", 
                                   fill = TRUE, 
                                   header = TRUE)


## Importamos los archivos csv de target
TargetAntimicrobiano <- read.table("./targets/Target_antimicrobiano.csv",
                                    sep = ",", 
                                    fill = TRUE, 
                                    header = TRUE)
TargetAntimicrobiano$Resultado <- str_replace(TargetAntimicrobiano$Resultado, "True", "TRUE")
TargetAntimicrobiano$Resultado <- str_replace(TargetAntimicrobiano$Resultado, "False", "FALSE")

# TargetAntimicrobiano$Resultado <- as.integer(as.logical(TargetAntimicrobiano$Resultado))
# TargetAntimicrobiano$comp <- ifelse(TargetAntimicrobiano$Molecula %in% DescriptoresLipinski$X, 1, NA)

## Limpieza de datos y tratamiento de los dataframes
Descriptores2D <- Descriptores2D %>%
  group_by(X) %>%
  filter(!duplicated(X))
Descriptores3D <- Descriptores3D %>%
  group_by(X) %>%
  filter(!duplicated(X))
DescriptoresLipinski <- DescriptoresLipinski %>%
  group_by(X) %>%
  filter(!duplicated(X))
TargetAntimicrobiano <- TargetAntimicrobiano %>%
  group_by(Molecula) %>%
  filter(!duplicated(Molecula))
Descriptores_Dataframe <- merge(Descriptores2D, 
                                Descriptores3D, 
                                by = "X", 
                                all = TRUE)
Descriptores_Dataframe <- merge(Descriptores_Dataframe, 
                                DescriptoresLipinski, 
                                by = "X", 
                                all = TRUE)
Descriptores_Dataframe_antimicrobianos <- merge(Descriptores_Dataframe,
                                                 TargetAntimicrobiano,
                                                 by.x = "X", 
                                                 all.x = TRUE, 
                                                 by.y = "Molecula", 
                                                 all.y = TRUE)

rownames(Descriptores_Dataframe_antimicrobianos) <- Descriptores_Dataframe_antimicrobianos$X
Descriptores_Dataframe_antimicrobianos <- Descriptores_Dataframe_antimicrobianos[!grepl("X", names(Descriptores_Dataframe_antimicrobianos))]


dim(Descriptores_Dataframe_antimicrobianos)
glimpse(Descriptores_Dataframe_antimicrobianos)
table(Descriptores_Dataframe_antimicrobianos$Resultado)
prop.table(table(Descriptores_Dataframe_antimicrobianos$Resultado)) %>% round(digits = 2)

# NIVEL BASAL 0.72

### Pre procesamiento
# Algunas no son utiles depndiendo el modelo solo para algunos modelos

nzv <- nearZeroVar(Descriptores_Dataframe_antimicrobianos)
datos_nzv <- Descriptores_Dataframe_antimicrobianos %>% nearZeroVar(saveMetrics = TRUE)
datos_nzv
Descriptores_df_am_filtrado <- Descriptores_Dataframe_antimicrobianos[, -nzv]

df_am_cor <- cor(Descriptores_df_am_filtrado[, -which(names(Descriptores_df_am_filtrado) == "Resultado")])
df_am_highlyCor <- findCorrelation(df_am_cor, cutoff = .75)
Descriptores_df_am_filtrado <- Descriptores_df_am_filtrado[, -df_am_highlyCor]

# va a funcionar en este caso si no se aplica el de highly cor, al menos para esta base de datos

#comboInfo <- findLinearCombos(Descriptores_df_am_filtrado[, -which(names(Descriptores_df_am_filtrado) == "Resultado")])
#Descriptores_df_am_filtrado <- Descriptores_df_am_filtrado[, -comboInfo$remove]

# Estableciendo reproducibilidad
set.seed(54)

## Data splitting
train_index <- createDataPartition(y = Descriptores_df_am_filtrado$Resultado,
                                   p = 0.75, 
                                   list = FALSE)

df_am_entrenamiento <- Descriptores_df_am_filtrado[train_index, ]
df_am_test <- Descriptores_df_am_filtrado[-train_index, ]

prop.table(table(df_am_entrenamiento$Resultado))

prop.table(table(df_am_test$Resultado))

# Mas metodos de preprocesado
objeto_recipe <- recipe(formula = Resultado ~ ., data = df_am_entrenamiento)
objeto_recipe
objeto_recipe <- objeto_recipe %>%  step_center(all_numeric())
objeto_recipe <- objeto_recipe %>%  step_scale(all_numeric())
trained_recipe <- prep(objeto_recipe, training = df_am_entrenamiento)
trained_recipe

df_am_entrenamiento <- bake(trained_recipe, new_data = df_am_entrenamiento)
df_am_test <- bake(trained_recipe, new_data = df_am_test)

glimpse(df_am_entrenamiento)

### Seleccion de predictores
## Eliminacion recursiva de variables (rfe)
subsets <- (1:148)
repeticiones <- 30
seeds <- vector(mode = "list", length = repeticiones + 1)
for (i in 1:repeticiones) {
  seeds[[i]] <- sample.int(1000, length(subsets))
}
seeds[[repeticiones + 1]] <- sample.int(1000, 1)

ctrl_rfe_am <- rfeControl(functions = rfFuncs,
                       method = "boot",
                       number = repeticiones,
                       returnResamp = "all",
                       allowParallel = TRUE,
                       verbose = FALSE,
                       seeds = seeds)

rf_rfe_am <- rfe(Resultado ~ .,
              data = df_am_entrenamiento,
              sizes = subsets,
              metric = "Accuracy",
              rfeControl = ctrl_rfe_am,
              ntree = 500)

rf_rfe_am

rf_rfe_am$optVariables
rf_rfe_am$resample %>% head(30)

rf_rfe_am$resample %>% group_by(Variables) %>%
                    summarise(media_accuracy = mean(Accuracy),
                              media_kappa = mean(Kappa)) %>%
                    arrange(desc(media_accuracy))

ggplot(data = rf_rfe_am$results, aes(x = Variables, y = Accuracy)) +
  geom_line() +
  scale_x_continuous(breaks  = unique(rf_rfe_am$results$Variables)) +
  geom_point() +
  geom_errorbar(aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD),
                width = 0.2) +
  geom_point(data = rf_rfe_am$results %>% slice(which.max(Accuracy)),
             color = "green") +
  theme(axis.text.x= element_text(size = 7, angle = 75)) + 
  labs(y = "precisión") 
  theme_bw()

head(rf_rfe_am$variables, 10)

rf_rfe_am$variables %>% filter(Variables == 124) %>% group_by(var) %>%
                     summarise(influencia_promedio = mean(Overall),
                               influencia_sd = sd(Overall)) %>%
                     arrange(desc(influencia_promedio))

rf_rfe_am$variables %>% filter(Variables == 147) %>% group_by(var) %>%
                     summarise(influencia_promedio = mean(Overall),
                               influencia_sd = sd(Overall)) %>%
                     arrange(desc(influencia_promedio))
## Algoritmo genético
ga_ctrl_am <- gafsControl(functions = rfGA,
                       method = "repeatedcv",
                       repeats = 5,
                       allowParallel = TRUE,
                       genParallel = TRUE,
                       verbose = FALSE)

rf_ga_am <- gafs(x =  df_am_entrenamiento[, !names(df_am_entrenamiento) %in% c("Resultado")],
              y = df_am_entrenamiento$Resultado,
              iters = 20,
              popSize = 100,
              gafsControl = ga_ctrl_am,
              ntree = 100)
rf_ga_am
rf_ga_am$optVariables
rf_ga_am$external %>% group_by(Iter) %>% summarize(Precision_media = mean(Accuracy))

plot(rf_ga_am) +
  labs(x = "Generación",y = "Media", color = "Estimado") + 
  scale_color_manual(labels = c("Externo", "Interno"), values = c("red", "blue")) +
  theme_bw()

poblacion_generada <- function(popSize, n_variables, n_max){
  poblacion <- matrix(data = NA, nrow = popSize, ncol = n_variables)
  for(i in 1:popSize){
    numero_1s <- sample(x = 1:n_max, size = 1)
    individuo <- rep(0, times = n_variables)
    individuo[sample(x = 1:n_variables, size = numero_1s)] <- 1
    poblacion[i,] <- individuo
  }
  return(poblacion)
}

poblacion_inicial <- poblacion_generada(popSize = 100,
                                        n_variables = ncol(df_am_entrenamiento) - 1 ,
                                        n_max = 4)
poblacion_inicial

rf_ga2_am <- gafs(x = df_am_entrenamiento %>% select(-Resultado),
              y = df_am_entrenamiento$Resultado,
              iters = 20,
              popSize = 100,
              suggestions = poblacion_inicial,
              gafsControl = ga_ctrl_am,
              ntree = 100)
rf_ga2_am

plot(rf_ga2_am) +
  labs(x = "Generación",y = "Media", color = "Estimado") + 
  scale_color_manual(labels = c("Externo", "Interno"), values = c("red", "blue")) +
  theme_bw()

## Método de Filtrado

particiones <- 10
repeticiones <- 5
seeds <- sample.int(1000, particiones * repeticiones + 1)
ctrl_filtrado_am <- sbfControl(functions = rfSBF,
                            method = "repeatedcv",
                            number = 10,
                            repeats = 5,
                            seeds = seeds,
                            verbose = FALSE,
                            saveDetails = TRUE,
                            allowParallel = TRUE)

rf_sbf <- sbf(Resultado ~ .,
              data = df_am_entrenamiento,
              sbfControl = ctrl_filtrado_am,
              ntree = 1000)

rf_sbf

rf_sbf$optVariables

# metodo de recocido simulado 

safs_ctrl_am = safsControl(functions = caretSA,
                        method = "repeatedcv",
                        repeats = 5,
                        improve =50,
                        allowParallel = TRUE,
                        verbose = FALSE)

rf_safs_am <- safs(x = df_am_entrenamiento %>% select(-Resultado),
                y = df_am_entrenamiento$Resultado,
                iters = 10,
                safsControl = safs_ctrl_am)

rf_safs_am

plot(rf_safs_am) +
  labs(x = "Iteración",y = "Media", color = "Estimado") + 
  scale_color_manual(labels = c("Externo", "Interno"), values = c("red", "blue")) +
  theme_bw()


# Metodo seleccionado de feature selection
predictores_filtrados <- rf_ga_am$optVariables
df_am_entrenamiento <- df_am_entrenamiento[, c(predictores_filtrados, 'Resultado')]
df_am_entrenamiento
df_am_test <- df_am_test[, c(predictores_filtrados, 'Resultado')]
df_am_test

### Entrenamiento de modelos
set.seed(54)

particiones  <- 10
repeticiones <- 5

cv_index <- createFolds(y = df_am_entrenamiento$Resultado,
                        k = 20,
                        list = TRUE, 
                        returnTrain = TRUE)

fitControl <- trainControl(index = cv_index,
                           method = "repeatedcv",
                           number = particiones,
                           repeats = repeticiones,
                           returnResamp = "all",
                           search = "random",
                           verboseIter = FALSE,
                           allowParallel = TRUE,
                           savePredictions = 'final')

fitControl_adaptive <- trainControl(method = "adaptive_cv",
                                       number = particiones,
                                       repeats = repeticiones,
                                        adaptive = list(min = 5, alpha = 0.05,
                                                         method = "gls", complete = TRUE),
                                       returnResamp = "all",
                                       verboseIter = FALSE,
                                       allowParallel = TRUE,
                                       savePredictions = 'final')

fitControl_final <- trainControl(method = "adaptive_cv",
                                    number = particiones,
                                    repeats = repeticiones,
                                    adaptive = list(min = 5, alpha = 0.05,
                                                    method = "gls", complete = TRUE),
                                    returnResamp = "final",
                                    verboseIter = FALSE,
                                    allowParallel = TRUE,
                                    savePredictions = 'final')

fitControl_final_random <- trainControl(index = cv_index,
                                        method = "repeatedcv",
                                        number = particiones,
                                        repeats = repeticiones,
                                        returnResamp = "final",
                                        search = "random",
                                        verboseIter = FALSE,
                                        allowParallel = TRUE,
                                        savePredictions = 'final')

ROC_control <- trainControl(method = "repeatedcv",
                            number = particiones,
                            repeats = repeticiones,
                            returnResamp = "final",
                            verboseIter = FALSE,
                            summaryFunction = twoClassSummary,
                            classProbs = TRUE,
                            allowParallel = TRUE)


## "L2 Regularized Linear Support Vector Machines with Class Weights")
# grid random
set.seed(54)

svmLinearWeights2_df_am_final <- caret::train(Resultado ~ .,
                                             data = df_am_entrenamiento,
                                             method = "svmLinearWeights2",
                                             tuneLength = 25,
                                             trControl = fitControl_final)

svmLinearWeights2_df_am_final$finalModel
svmLinearWeights2_df_am_final
svmLinearWeights2_df_am_final$resample %>% head(10)
summary(svmLinearWeights2_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmLinearWeights2 <- ggplot(svmLinearWeights2_df_am_final,  highlight = TRUE) +
  labs(title = "svmLinearWeights2") +
  theme_bw()

grafica_hiperparametros_svmLinearWeights2

svmLinearWeights2_p1 <- ggplot(data = svmLinearWeights2_df_am_final$resample, aes(x = Accuracy)) +
      geom_density(alpha = 0.5, fill = "green") +
      geom_vline(xintercept = mean(svmLinearWeights2_df_am_final$resample$Accuracy),
                 linetype = "dashed") +
      labs(x = "precisión", y = "densidad")
      theme_bw()

svmLinearWeights2_p2  <- ggplot(data = svmLinearWeights2_df_am_final$resample, aes(x = 1, y = Accuracy)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
      geom_jitter(width = 0.05) +
      labs(y = "precisión") +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmLinearWeights2 <- ggarrange(svmLinearWeights2_p1, svmLinearWeights2_p2)
final_plot_svmLinearWeights2<- annotate_figure(
                final_plot_svmLinearWeights2,
                top = text_grob("svmLinearWeights2", size = 15))
final_plot_svmLinearWeights2

# No funciona ROC en este modelo
# ROC_svmLinearWeights2<- train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "svmLinearWeights2",
#                               metric = "ROC",
#                               trControl = ROC_control,
#                               family = "binomial")

# ROC_svmLinearWeights2

# predicciones_roc_svmLinearWeights2 <- predict(object = ROC_svmLinearWeights2,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_svmLinearWeights2 <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_svmLinearWeights2$TRUE)

# plot(curva_ROC_svmLinearWeights2)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_svmLinearWeights2, conf.level = 0.95)

#predicciones_probabilidad_svmLinearWeights2 <- predict(svmLinearWeights2_df_am_final, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_svmLinearWeights2 <- predict(svmLinearWeights2_df_am_final,
                                                 newdata = df_am_test,
                                                 type = "raw")
predicciones_crudas_svmLinearWeights2 %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmLinearWeights2, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmLinearWeights2, df_am_test$Resultado,
           expected=TRUE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,   
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmLinearWeights2 != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Linear Support Vector Machines with Class Weights")
# grid random
#problemas de potencia returnresamp final

svmLinearWeights_df_am_final <- train(Resultado ~ .,
                                       data = df_am_entrenamiento,
                                       method = "svmLinearWeights",
                                       tuneLength = 15,
                                       trControl = fitControl_final)

svmLinearWeights_df_am_final$finalModel
svmLinearWeights_df_am_final
svmLinearWeights_df_am_final$resample %>% head(10)
summary(svmLinearWeights_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmLinearWeights <- ggplot(svmLinearWeights_df_am_final,  highlight = TRUE) +
  labs(title = "svmLinearWeights") +
  theme_bw()

svmLinearWeights_p1 <- ggplot(data = svmLinearWeights_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(svmLinearWeights_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmLinearWeights_p2  <- ggplot(data = svmLinearWeights_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmLinearWeights <- ggarrange(svmLinearWeights_p1, svmLinearWeights_p2)
final_plot_svmLinearWeights <- annotate_figure(
  final_plot_svmLinearWeights,
  top = text_grob("svmLinearWeights", size = 15))
final_plot_svmLinearWeights

ROC_svmLinearWeights <- train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "svmLinearWeights",
                              trControl = ROC_control,
                              metric = "ROC",
                             family = "binomial")

ROC_svmLinearWeights

predicciones_roc_svmLinearWeights <- predict(object = ROC_svmLinearWeights,
                           newdata = df_am_test,
                           type = "prob")
predicciones_roc_svmLinearWeights
curva_ROC_svmLinearWeights <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_svmLinearWeights$TRUE.)

grafica_ROC_svmLinearWeights <- plot(curva_ROC_svmLinearWeights, main = "svmLinearWeights") 

# Intervalo de confianza de la curva
ci.auc(curva_ROC_svmLinearWeights, conf.level = 0.95)

predicciones_crudas_svmLinearWeights <- predict(svmLinearWeights_df_am_final,
                                                 newdata = df_am_test,
                                                 type = "raw")
predicciones_crudas_svmLinearWeights %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmLinearWeights, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmLinearWeights, df_am_test$Resultado,
            expected=FALSE, prop.r=TRUE, prop.c=TRUE,
            prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
            resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")


error_test <- mean(predicciones_crudas_svmLinearWeights != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Support Vector Machines with Class Weights")
# adaptive

svmRadialWeights_df_am_final <- train(Resultado ~ .,
                                      data = df_am_entrenamiento,
                                      method = "svmRadialWeights",
                                      tuneLength = 40,
                                      trControl = fitControl_final)

svmRadialWeights_df_am_final$finalModel
svmRadialWeights_df_am_final
svmRadialWeights_df_am_final$resample %>% head(10)
summary(svmRadialWeights_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmRadialWeights <- ggplot(svmRadialWeights_df_am_final,  highlight = TRUE) +
  labs(title = "svmRadialWeights") +
  theme_bw()

svmRadialWeights_p1 <- ggplot(data = svmRadialWeights_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(svmRadialWeights_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmRadialWeights_p2  <- ggplot(data = svmRadialWeights_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmRadialWeights <- ggarrange(svmRadialWeights_p1, svmRadialWeights_p2)
final_plot_svmRadialWeights <- annotate_figure(
  final_plot_svmRadialWeights,
  top = text_grob("svmRadialWeights", size = 15))
final_plot_svmRadialWeights

ROC_svmRadialWeights <- train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "svmRadialWeights",
                              trControl = ROC_control,
                              metric = "ROC",
                             family = "binomial")

ROC_svmRadialWeights

predicciones_roc_svmRadialWeights <- predict(object = ROC_svmRadialWeights,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_svmRadialWeights <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_svmRadialWeights$TRUE.)

grafica_ROC_svmRadialWeights <- plot(curva_ROC_svmRadialWeights, main = "svmRadialWeights")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_svmRadialWeights, conf.level = 0.95)

predicciones_crudas_svmRadialWeights <- predict(svmRadialWeights_df_am_final,
                                                newdata = df_am_test,
                                                type = "raw")
predicciones_crudas_svmRadialWeights %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmRadialWeights, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmRadialWeights, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmRadialWeights != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "L2 Regularized Support Vector Machine (dual) with Linear Kernel")
# Random

svmLinear3_df_am_final <- train(Resultado ~ .,
                                      data = df_am_entrenamiento,
                                      method = "svmLinear3",
                                      tuneLength = 40,
                                      trControl = fitControl_final)

svmLinear3_df_am_final$finalModel
svmLinear3_df_am_final
svmLinear3_df_am_final$resample %>% head(10)
summary(svmLinear3_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmLinear3 <- ggplot(svmLinear3_df_am_final,  highlight = TRUE) +
  labs(title = "svmLinear3") +
  theme_bw()

svmLinear3_p1 <- ggplot(data = svmLinear3_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(svmLinear3_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmLinear3_p2  <- ggplot(data = svmLinear3_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmLinear3 <- ggarrange(svmLinear3_p1, svmLinear3_p2)
final_plot_svmLinear3 <- annotate_figure(
  final_plot_svmLinear3,
  top = text_grob("svmLinear3", size = 15))
final_plot_svmLinear3

# No funciona ROC
# ROC_svmLinear3 <- train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "svmLinear3",
#                               trControl = ROC_control,
#                               metric = "ROC",
#                              family = "binomial")

# ROC_svmLinear3

# predicciones_roc_svmLinear3 <- predict(object = ROC_svmLinear3,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_svmLinear3 <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_svmLinear3$TRUE)

# plot(curva_ROC_svmLinear3)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_svmLinear3, conf.level = 0.95)

# predicciones_probabilidad_svmLinear3 <- predict(ROC_svmLinear3, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_svmLinear3 <- predict(svmLinear3_df_am_final,
                                                newdata = df_am_test,
                                                type = "raw")
predicciones_crudas_svmLinear3 %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmLinear3, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmLinear3, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmLinear3 != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Least Squares Support Vector Machine with Polynomial Kernel")
# Da valores raros con NaNs
# adaptive
#problemas de potencia returnresamp final
# 
# lssvmPoly_final <- train(Resultado ~ .,
#                                 data = df_am_entrenamiento,
#                                 method = "lssvmPoly",
#                                 tuneLength = 30,
#                                 trControl = fitControl_final)
# 
# lssvmPoly_final$finalModel
# lssvmPoly_final
# lssvmPoly_final$resample %>% head(10)
# summary(lssvmPoly_final$resample$Accuracy)
# 
# grafica_hiperparametros_lssvmPoly <- ggplot(lssvmPoly_final,  highlight = TRUE) +
#   labs(title = "lssvmPoly") +
#   theme_bw()
# 
# lssvmPoly_p1 <- ggplot(data = lssvmPoly_final$resample, aes(x = Accuracy)) +
#   geom_density(alpha = 0.5, fill = "green") +
#   geom_vline(xintercept = mean(lssvmPoly_final$resample$Accuracy),
#              linetype = "dashed") +
#   labs(x = "precisión", y = "densidad")
# theme_bw()
# 
# lssvmPoly_p2  <- ggplot(data = lssvmPoly_final$resample, aes(x = 1, y = Accuracy)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
#   geom_jitter(width = 0.05) +
#   labs(y = "precisión") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# final_plot_lssvmPoly <- ggarrange(lssvmPoly_p1, lssvmPoly_p2)
# final_plot_lssvmPoly <- annotate_figure(
#   final_plot_lssvmPoly,
#   top = text_grob("lssvmPoly", size = 15))
# final_plot_lssvmPoly
# 
# # No funciona ROC
# # ROC_lssvmPoly <- train(make.names(Resultado) ~ .,
# #                               data = df_am_entrenamiento,
# #                               method = "lssvmPoly",
# #                               trControl = ROC_control,
# #                               metric = "ROC",
# #                              family = "binomial")
# 
# # ROC_lssvmPoly
# 
# # predicciones_roc_lssvmPoly <- predict(object = ROC_lssvmPoly,
# #                            newdata = df_am_test,
# #                            type = "prob")
# # curva_ROC_lssvmPoly <- roc(response = df_am_test$Resultado,
# #                                      predictor = predicciones_roc_lssvmPoly$TRUE.)
# 
# # plot(curva_ROC_lssvmPoly)
# 
# # Intervalo de confianza de la curva
# # ci.auc(curva_ROC_lssvmPoly, conf.level = 0.95)
# 
# # predicciones_probabilidad_lssvmPoly <- predict(lssvmPoly_final, newdata = df_am_test,
# #                                               type = "prob")
# 
# predicciones_crudas_lssvmPoly <- predict(lssvmPoly_final,
#                                           newdata = df_am_test,
#                                           type = "raw")
# predicciones_crudas_lssvmPoly %>% head(10)
# 
# caret::confusionMatrix(data = predicciones_crudas_lssvmPoly, reference = df_am_test$Resultado,
#                 positive = "TRUE")
# 
# CrossTable(predicciones_crudas_lssvmPoly, df_am_test$Resultado,
#            expected=FALSE, prop.r=TRUE, prop.c=TRUE,
#            prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
#            resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")
# 
# error_test <- mean(predicciones_crudas_lssvmPoly != df_am_test$Resultado)
# paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Least Squares Support Vector Machine with Radial Basis Function Kernel")
# Adaptive

lssvmRadial_df_am_final <- train(Resultado ~ .,
                           data = df_am_entrenamiento,
                           method = "lssvmRadial",
                           trControl = fitControl_final,
                           tuneLenght = 40)

lssvmRadial_df_am_final$finalModel
lssvmRadial_df_am_final
lssvmRadial_df_am_final$resample %>% head(10)
summary(lssvmRadial_df_am_final$resample$Accuracy)

grafica_hiperparametros_lssvmRadial <- ggplot(lssvmRadial_df_am_final,  highlight = TRUE) +
  labs(title = "lssvmRadial") +
  theme_bw()

lssvmRadial_p1 <- ggplot(data = lssvmRadial_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(lssvmRadial_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

lssvmRadial_p2  <- ggplot(data = lssvmRadial_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_lssvmRadial <- ggarrange(lssvmRadial_p1, lssvmRadial_p2)
final_plot_lssvmRadial <- annotate_figure(
  final_plot_lssvmRadial,
  top = text_grob("lssvmRadial", size = 15))
final_plot_lssvmRadial

# No funciona ROC
# ROC_lssvmRadial <- train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "lssvmRadial",
#                               trControl = ROC_control,
#                               metric = "ROC",
#                              family = "binomial")

# ROC_lssvmRadial

# predicciones_roc_lssvmRadial <- predict(object = ROC_lssvmRadial,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_lssvmRadial <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_lssvmRadial$TRUE)

# plot(curva_ROC_lssvmRadial)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_lssvmRadial, conf.level = 0.95)

# predicciones_probabilidad_lssvmRadial <- predict(lssvmRadial_df_am_final, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_lssvmRadial <- predict(lssvmRadial_df_am_final,
                                         newdata = df_am_test,
                                         type = "raw")
predicciones_crudas_lssvmRadial %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_lssvmRadial, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_lssvmRadial, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_lssvmRadial != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")
# 
# ## "Support Vector Machines with Linear Kernel
# Los hiperparametros no son aplicables al adaptive cv.
# Da un numero desigual de resamples al utilizar repeated cv, por lo que dificulta su comparacion
# El resultado no es muy bueno un, 67%
# 
# fitControl_svmLinear_final <- trainControl(index = cv_index,
#                                      method = "repeatedcv",
#                                      number = particiones,
#                                      repeats = repeticiones,
#                                      returnResamp = "final",
#                                      verboseIter = FALSE,
#                                      allowParallel = TRUE,
#                                      savePredictions = 'final')
# 
# svmLinear_df_am_final <- train(Resultado ~ .,
#                          data = df_am_entrenamiento,
#                          method = "svmLinear",
#                          trControl = fitControl_svmLinear_final,
#                          tuneLenght = 40)
# 
# svmLinear_df_am_final$finalModel
# svmLinear_df_am_final
# svmLinear_df_am_final$resample %>% head(10)
# summary(svmLinear_df_am_final$resample$Accuracy)
# 
# # no tiene muchos parametros como para generar una grafica
# # grafica_hiperparametros_svmLinear <- ggplot(svmLinear_df_am_final,  highlight = TRUE) +
# #   labs(title = "svmLinear") +
# #   theme_bw()
# 
# svmLinear_p1 <- ggplot(data = svmLinear_df_am_final$resample, aes(x = Accuracy)) +
#   geom_density(alpha = 0.5, fill = "green") +
#   geom_vline(xintercept = mean(svmLinear_df_am_final$resample$Accuracy),
#              linetype = "dashed") +
#   labs(x = "precisión", y = "densidad")
# theme_bw()
# 
# svmLinear_p2  <- ggplot(data = svmLinear_df_am_final$resample, aes(x = 1, y = Accuracy)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
#   geom_jitter(width = 0.05) +
#   labs(y = "precisión") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# final_plot_svmLinear <- ggarrange(svmLinear_p1, svmLinear_p2)
# final_plot_svmLinear <- annotate_figure(
#   final_plot_svmLinear,
#   top = text_grob("svmLinear", size = 15))
# final_plot_svmLinear
# 
# ROC_svmLinear <- train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "svmLinear",
#                               trControl = ROC_control,
#                               metric = "ROC",
#                              family = "binomial")
# 
# ROC_svmLinear
# 
# predicciones_roc_svmLinear <- predict(object = ROC_svmLinear,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_svmLinear <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_svmLinear$TRUE.)
# 
# grafica_ROC_svmLinear <- plot(curva_ROC_svmLinear, main = "svmLinear")
# 
# # Intervalo de confianza de la curva
# ci.auc(curva_ROC_svmLinear, conf.level = 0.95)
# 
# 
# predicciones_crudas_svmLinear <- predict(svmLinear_df_am_final,
#                                            newdata = df_am_test,
#                                            type = "raw")
# predicciones_crudas_svmLinear %>% head(10)
# 
# caret::confusionMatrix(data = predicciones_crudas_svmLinear, reference = df_am_test$Resultado,
#                 positive = "TRUE")
# 
# CrossTable(predicciones_crudas_svmLinear, df_am_test$Resultado,
#            expected=FALSE, prop.r=TRUE, prop.c=TRUE,
#            prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
#            resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")
# 
# error_test <- mean(predicciones_crudas_svmLinear != df_am_test$Resultado)
# paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Support Vector Machines with Linear Kernel")
# adaptive

svmLinear2_df_am_final <- train(Resultado ~ .,
                          data = df_am_entrenamiento,
                          method = "svmLinear2",
                          trControl = fitControl_final,
                          tuneLenght = 40)

svmLinear2_df_am_final$finalModel
svmLinear2_df_am_final
svmLinear2_df_am_final$resample %>% head(10)
summary(svmLinear2_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmLinear2 <- ggplot(svmLinear2_df_am_final,  highlight = TRUE) +
  labs(title = "svmLinear2") +
  theme_bw()

svmLinear2_p1 <- ggplot(data = svmLinear2_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(svmLinear2_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmLinear2_p2  <- ggplot(data = svmLinear2_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmLinear2 <- ggarrange(svmLinear2_p1, svmLinear2_p2)
final_plot_svmLinear2 <- annotate_figure(
  final_plot_svmLinear2,
  top = text_grob("svmLinear2", size = 15))
final_plot_svmLinear2

ROC_svmLinear2 <- train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "svmLinear2",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_svmLinear2

predicciones_roc_svmLinear2 <- predict(object = ROC_svmLinear2,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_svmLinear2 <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_svmLinear2$TRUE.)

grafica_ROC_svmLinear2 <- plot(curva_ROC_svmLinear2, main = "svmLinear2")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_svmLinear2, conf.level = 0.95)


predicciones_crudas_svmLinear2 <- predict(svmLinear2_df_am_final,
                                         newdata = df_am_test,
                                         type = "raw")
predicciones_crudas_svmLinear2 %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmLinear2, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmLinear2, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmLinear2 != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Support Vector Machines with Polynomial kernel")
#adaptive

svmPoly_df_am_final <- train(Resultado ~ .,
                       data = df_am_entrenamiento,
                       method = "svmPoly",
                       trControl = fitControl_final,
                       tuneLenght = 40)

svmPoly_df_am_final$finalModel
svmPoly_df_am_final
svmPoly_df_am_final$resample %>% head(10)
summary(svmPoly_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmPoly <- ggplot(svmPoly_df_am_final,  highlight = TRUE) +
  labs(title = "svmPoly") +
  theme_bw()

svmPoly_p1 <- ggplot(data = svmPoly_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept =mean(svmPoly_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmPoly_p2  <- ggplot(data = svmPoly_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmPoly <- ggarrange(svmPoly_p1, svmPoly_p2)
final_plot_svmPoly <- annotate_figure(
  final_plot_svmPoly,
  top = text_grob("svmPoly", size = 15))
final_plot_svmPoly

ROC_svmPoly<- train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "svmPoly",
                              trControl = ROC_control,
                              metric = "ROC",
                             family = "binomial")

ROC_svmPoly

predicciones_roc_svmPoly <- predict(object = ROC_svmPoly,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_svmPoly <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_svmPoly$TRUE.)

grafica_ROC_svmPoly <- plot(curva_ROC_svmPoly, main = "svmPoly")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_svmPoly, conf.level = 0.95)

predicciones_crudas_svmPoly <- predict(svmPoly_df_am_final,
                                          newdata = df_am_test,
                                          type = "raw")
predicciones_crudas_svmPoly %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmPoly, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmPoly, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmPoly != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Support Vector Machines with Radial Basis Function Kernel")
# adaptative

svmRadial_df_am_final <- train(Resultado ~ .,
                         data = df_am_entrenamiento,
                         method = "svmRadial",
                         trControl = fitControl_final,
                         tuneLenght = 40)

svmRadial_df_am_final$finalModel
svmRadial_df_am_final
svmRadial_df_am_final$resample %>% head(10)
summary(svmRadial_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmRadial <- ggplot(svmRadial_df_am_final,  highlight = TRUE) +
  labs(title = "svmRadial") +
  theme_bw()

svmRadial_p1 <- ggplot(data = svmRadial_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(svmRadial_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmRadial_p2  <- ggplot(data = svmRadial_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmRadial <- ggarrange(svmRadial_p1, svmRadial_p2)
final_plot_svmRadial <- annotate_figure(
  final_plot_svmRadial,
  top = text_grob("svmRadial", size = 15))
final_plot_svmRadial

ROC_svmRadial<- train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "svmRadial",
                              trControl = ROC_control,
                             family = "binomial", 
                              metric = "ROC")

ROC_svmRadial

predicciones_roc_svmRadial <- predict(object = ROC_svmRadial,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_svmRadial <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_svmRadial$TRUE.)

grafica_ROC_svmRadial <- plot(curva_ROC_svmRadial, main = "svmRadial")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_svmRadial, conf.level = 0.95)

predicciones_crudas_svmRadial <- predict(svmRadial_df_am_final,
                                       newdata = df_am_test,
                                       type = "raw")
predicciones_crudas_svmRadial %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmRadial, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmRadial, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmRadial != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Support Vector Machines with Radial Basis Function Kernel")
#Random

svmRadialCost_df_am_final <- train(Resultado ~ .,
                             data = df_am_entrenamiento,
                             method = "svmRadialCost",
                             trControl = fitControl_final,
                             tuneLenght = 40)

svmRadialCost_df_am_final$finalModel
svmRadialCost_df_am_final
svmRadialCost_df_am_final$resample %>% head(10)
summary(svmRadialCost_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmRadialCost <- ggplot(svmRadialCost_df_am_final,  highlight = TRUE) +
  labs(title = "svmRadialCost") +
  theme_bw()

svmRadialCost_p1 <- ggplot(data = svmRadialCost_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(svmRadialCost_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmRadialCost_p2  <- ggplot(data = svmRadialCost_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmRadialCost <- ggarrange(svmRadialCost_p1, svmRadialCost_p2)
final_plot_svmRadialCost <- annotate_figure(
  final_plot_svmRadialCost,
  top = text_grob("svmRadialCost", size = 15))
final_plot_svmRadialCost

ROC_svmRadialCost<- train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "svmRadialCost",
                              trControl = ROC_control,
                              metric = "ROC",
                             family = "binomial")

ROC_svmRadialCost

predicciones_roc_svmRadialCost <- predict(object = ROC_svmRadialCost,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_svmRadialCost <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_svmRadialCost$TRUE.)

grafica_ROC_svmRadialCost <- plot(curva_ROC_svmRadialCost, main = "svmRadialCost")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_svmRadialCost, conf.level = 0.95)

predicciones_crudas_svmRadialCost <- predict(svmRadialCost_df_am_final,
                                         newdata = df_am_test,
                                         type = "raw")
predicciones_crudas_svmRadialCost %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmRadialCost, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmRadialCost, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmRadialCost != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Support Vector Machines with Radial Basis Function Kernel")
# Random

svmRadialSigma_df_am_final <- train(Resultado ~ .,
                              data = df_am_entrenamiento,
                              method = "svmRadialSigma",
                              trControl = fitControl_final,
                              tuneLenght = 40)

svmRadialSigma_df_am_final$finalModel
svmRadialSigma_df_am_final
svmRadialSigma_df_am_final$resample %>% head(10)
summary(svmRadialSigma_df_am_final$resample$Accuracy)

grafica_hiperparametros_svmRadialSigma <- ggplot(svmRadialSigma_df_am_final,  highlight = TRUE) +
  labs(title = "svmRadialSigma") +
  theme_bw()

svmRadialSigma_p1 <- ggplot(data = svmRadialSigma_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(svmRadialSigma_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

svmRadialSigma_p2  <- ggplot(data = svmRadialSigma_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_svmRadialSigma <- ggarrange(svmRadialSigma_p1, svmRadialSigma_p2)
final_plot_svmRadialSigma <- annotate_figure(
  final_plot_svmRadialSigma,
  top = text_grob("svmRadialSigma", size = 15))
final_plot_svmRadialSigma

ROC_svmRadialSigma<- train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "svmRadialSigma",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_svmRadialSigma

predicciones_roc_svmRadialSigma <- predict(object = ROC_svmRadialSigma,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_svmRadialSigma <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_svmRadialSigma$TRUE.)

grafica_ROC_svmRadialSigma <- plot(curva_ROC_svmRadialSigma, main = "svmRadialSigma")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_svmRadialSigma, conf.level = 0.95)

predicciones_crudas_svmRadialSigma <- predict(svmRadialSigma_df_am_final,
                                             newdata = df_am_test,
                                             type = "raw")
predicciones_crudas_svmRadialSigma %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_svmRadialSigma, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_svmRadialSigma, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_svmRadialSigma != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Multi-Layer Perceptron")
# adaptive

mlp_df_am_final <- caret::train(Resultado ~ .,
                          data = df_am_entrenamiento,
                          method = "mlp",
                          trControl = fitControl_final,
                          tuneLenght = 40)

mlp_df_am_final$finalModel
mlp_df_am_final
mlp_df_am_final$resample %>% head(10)
summary(mlp_df_am_final$resample$Accuracy)

grafica_hiperparametros_mlp <- ggplot(mlp_df_am_final,  highlight = TRUE) +
  labs(title = "mlp") +
  theme_bw()

mlp_p1 <- ggplot(data = mlp_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(mlp_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

mlp_p2  <- ggplot(data = mlp_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_mlp <- ggarrange(mlp_p1, mlp_p2)
final_plot_mlp <- annotate_figure(
  final_plot_mlp,
  top = text_grob("mlp", size = 15))
final_plot_mlp

ROC_mlp <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "mlp",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_mlp

predicciones_roc_mlp <- predict(object = ROC_mlp,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_mlp <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_mlp$TRUE.)

grafica_ROC_mlp <- plot(curva_ROC_mlp, main = "mlp")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_mlp, conf.level = 0.95)

predicciones_crudas_mlp <- predict(mlp_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_mlp %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_mlp, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_mlp, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_mlp != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Multi-layer Perceptron 2")
# adaptive

mlpWeightDecay_df_am_final <- caret::train(Resultado ~ .,
                                    data = df_am_entrenamiento,
                                    method = "mlpWeightDecay",
                                    trControl = fitControl_final,
                                    tuneLenght = 40)

mlpWeightDecay_df_am_final$finalModel
mlpWeightDecay_df_am_final
mlpWeightDecay_df_am_final$resample %>% head(10)
summary(mlpWeightDecay_df_am_final$resample$Accuracy)

grafica_hiperparametros_mlpWeightDecay <- ggplot(mlpWeightDecay_df_am_final,  highlight = TRUE) +
  labs(title = "mlpWeightDecay") +
  theme_bw()

mlpWeightDecay_p1 <- ggplot(data = mlpWeightDecay_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(mlpWeightDecay_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

mlpWeightDecay_p2  <- ggplot(data = mlpWeightDecay_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_mlpWeightDecay <- ggarrange(mlpWeightDecay_p1, mlpWeightDecay_p2)
final_plot_mlpWeightDecay <- annotate_figure(
  final_plot_mlpWeightDecay,
  top = text_grob("mlpWeightDecay", size = 15))
final_plot_mlpWeightDecay

ROC_mlpWeightDecay <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "mlpWeightDecay",
                              trControl = ROC_control,
                              metric = "ROC",
                             family = "binomial")

ROC_mlpWeightDecay

predicciones_roc_mlpWeightDecay <- predict(object = ROC_mlpWeightDecay,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_mlpWeightDecay <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_mlpWeightDecay$TRUE.)

grafica_ROC_mlpWeightDecay <- plot(curva_ROC_mlpWeightDecay, main = "mlpWeightDecay")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_mlpWeightDecay, conf.level = 0.95)

predicciones_probabilidad_mlpWeightDecay <- predict(mlpWeightDecay_df_am_final, newdata = df_am_test,
                                              type = "prob")

predicciones_crudas_mlpWeightDecay <- predict(mlpWeightDecay_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_mlpWeightDecay %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_mlpWeightDecay, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_mlpWeightDecay, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_mlpWeightDecay != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Multi-Layer Perceptron, multiple layers")
#adaptive

mlpWeightDecayML_df_am_final <- caret::train(Resultado ~ .,
                                      data = df_am_entrenamiento,
                                      method = "mlpWeightDecayML",
                                      trControl = fitControl_final,
                                      tuneLenght = 40)

mlpWeightDecayML_df_am_final$finalModel
mlpWeightDecayML_df_am_final
mlpWeightDecayML_df_am_final$resample %>% head(10)
summary(mlpWeightDecayML_df_am_final$resample$Accuracy)

grafica_hiperparametros_mlpWeightDecayML <- ggplot(mlpWeightDecayML_df_am_final,  highlight = TRUE) +
  labs(title = "mlpWeightDecayML") +
  theme_bw()

mlpWeightDecayML_p1 <- ggplot(data = mlpWeightDecayML_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(mlpWeightDecayML_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

mlpWeightDecayML_p2  <- ggplot(data = mlpWeightDecayML_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_mlpWeightDecayML <- ggarrange(mlpWeightDecayML_p1, mlpWeightDecayML_p2)
final_plot_mlpWeightDecayML <- annotate_figure(
  final_plot_mlpWeightDecayML,
  top = text_grob("mlpWeightDecayML", size = 15))
final_plot_mlpWeightDecayML

ROC_mlpWeightDecayML <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "mlpWeightDecayML",
                              trControl = ROC_control,
                              metric = "ROC",
                             family = "binomial")

ROC_mlpWeightDecayML

predicciones_roc_mlpWeightDecayML <- predict(object = ROC_mlpWeightDecayML,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_mlpWeightDecayML <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_mlpWeightDecayML$TRUE.)

grafica_ROC_mlpWeightDecayML <- plot(curva_ROC_mlpWeightDecayML, main = "mlpWeightDecayML")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_mlpWeightDecayML, conf.level = 0.95)

predicciones_crudas_mlpWeightDecayML <- predict(mlpWeightDecayML_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_mlpWeightDecayML %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_mlpWeightDecayML, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_mlpWeightDecayML, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_mlpWeightDecayML != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Multi-Layer Perceptron, with multiple layers")
# random

mlpML_df_am_final<- caret::train(Resultado ~ .,
                           data = df_am_entrenamiento,
                           method = "mlpML",
                           trControl = fitControl_final,
                           tuneLenght = 40)

mlpML_df_am_final$finalModel
mlpML_df_am_final
mlpML_df_am_final$resample %>% head(10)
summary(mlpML_df_am_final$resample$Accuracy)

grafica_hiperparametros_mlpML <- ggplot(mlpML_df_am_final,  highlight = TRUE) +
  labs(title = "mlpML") +
  theme_bw()

mlpML_p1 <- ggplot(data = mlpML_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(mlpML_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

mlpML_p2  <- ggplot(data = mlpML_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_mlpML <- ggarrange(mlpML_p1, mlpML_p2)
final_plot_mlpML <- annotate_figure(
  final_plot_mlpML,
  top = text_grob("mlpML", size = 15))
final_plot_mlpML

ROC_mlpML <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "mlpML",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_mlpML

predicciones_roc_mlpML <- predict(object = ROC_mlpML,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_mlpML <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_mlpML$TRUE.)

grafica_ROC_mlpML <- plot(curva_ROC_mlpML, main = "mlpML")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_mlpML, conf.level = 0.95)

predicciones_crudas_mlpML <- predict(mlpML_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_mlpML %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_mlpML, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_mlpML, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_mlpML != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Oblique Random Forest

ORFlog_df_am_final <- caret::train(Resultado ~ .,
                            data = df_am_entrenamiento,
                            method = "ORFlog",
                            trControl = fitControl_final,
                            tuneLenght = 40)

ORFlog_df_am_final$finalModel
ORFlog_df_am_final
ORFlog_df_am_final$resample %>% head(10)
summary(ORFlog_df_am_final$resample$Accuracy)

grafica_hiperparametros_ORFlog <- ggplot(ORFlog_df_am_final,  highlight = TRUE) +
  labs(title = "ORFlog") +
  theme_bw()

ORFlog_p1 <- ggplot(data = ORFlog_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(ORFlog_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

ORFlog_p2  <- ggplot(data = ORFlog_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_ORFlog <- ggarrange(ORFlog_p1, ORFlog_p2)
final_plot_ORFlog <- annotate_figure(
  final_plot_ORFlog,
  top = text_grob("ORFlog", size = 15))
final_plot_ORFlog

ROC_ORFlog <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "ORFlog",
                           metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_ORFlog

predicciones_roc_ORFlog <- predict(object = ROC_ORFlog,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_ORFlog <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_ORFlog$TRUE.)

grafica_ROC_ORFlog <- plot(curva_ROC_ORFlog, main = "ORFlog")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_ORFlog, conf.level = 0.95)

predicciones_crudas_ORFlog <- predict(ORFlog_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")

predicciones_crudas_ORFlog %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_ORFlog, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_ORFlog, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_ORFlog != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Oblique Random Forest

ORFpls_df_am_final <- caret::train(Resultado ~ .,
                            data = df_am_entrenamiento,
                            method = "ORFpls",
                            trControl = fitControl_final,
                            tuneLenght = 40)

ORFpls_df_am_final$finalModel
ORFpls_df_am_final
ORFpls_df_am_final$resample %>% head(10)
summary(ORFpls_df_am_final$resample$Accuracy)

grafica_hiperparametros_ORFpls <- ggplot(ORFpls_df_am_final,  highlight = TRUE) +
  labs(title = "ORFpls") +
  theme_bw()

ORFpls_p1 <- ggplot(data = ORFpls_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(ORFpls_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

ORFpls_p2  <- ggplot(data = ORFpls_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_ORFpls <- ggarrange(ORFpls_p1, ORFpls_p2)
final_plot_ORFpls <- annotate_figure(
  final_plot_ORFpls,
  top = text_grob("ORFpls", size = 15))
final_plot_ORFpls

ROC_ORFpls <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "ORFpls",
                           metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_ORFpls

predicciones_roc_ORFpls <- predict(object = ROC_ORFpls,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_ORFpls <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_ORFpls$TRUE.)

grafica_ROC_ORFpls <- plot(curva_ROC_ORFpls, main = "ORFpls")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_ORFpls, conf.level = 0.95)

predicciones_crudas_ORFpls <- predict(ORFpls_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_ORFpls %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_ORFpls, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_ORFpls, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_ORFpls != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Oblique Random Forest

ORFridge_df_am_final <- caret::train(Resultado ~ .,
                              data = df_am_entrenamiento,
                              method = "ORFridge",
                              trControl = fitControl_final,
                              tuneLenght = 40)

ORFridge_df_am_final$finalModel
ORFridge_df_am_final
ORFridge_df_am_final$resample %>% head(10)
summary(ORFridge_df_am_final$resample$Accuracy)

grafica_hiperparametros_ORFridge <- ggplot(ORFridge_df_am_final,  highlight = TRUE) +
  labs(title = "ORFridge") +
  theme_bw()

ORFridge_p1 <- ggplot(data = ORFridge_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(ORFridge_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

ORFridge_p2  <- ggplot(data = ORFridge_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_ORFridge <- ggarrange(ORFridge_p1, ORFridge_p2)
final_plot_ORFridge <- annotate_figure(
  final_plot_ORFridge,
  top = text_grob("ORFridge", size = 15))
final_plot_ORFridge

ROC_ORFridge <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "ORFridge",
                             metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_ORFridge

predicciones_roc_ORFridge <- predict(object = ROC_ORFridge,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_ORFridge <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_ORFridge$TRUE.)

grafica_ROC_ORFridge <- plot(curva_ROC_ORFridge, main = "ORFridge")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_ORFridge, conf.level = 0.95)

predicciones_crudas_ORFridge <- predict(ORFridge_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_ORFridge %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_ORFridge, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_ORFridge, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_ORFridge != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Oblique Random Forest

ORFsvm_df_am_final <- caret::train(Resultado ~ .,
                            data = df_am_entrenamiento,
                            method = "ORFsvm",
                            trControl = fitControl_final,
                            tuneLenght = 40)

ORFsvm_df_am_final$finalModel
ORFsvm_df_am_final
ORFsvm_df_am_final$resample %>% head(10)
summary(ORFsvm_df_am_final$resample$Accuracy)

grafica_hiperparametros_ORFsvm <- ggplot(ORFsvm_df_am_final,  highlight = TRUE) +
  labs(title = "ORFsvm") +
  theme_bw()

ORFsvm_p1 <- ggplot(data = ORFsvm_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(ORFsvm_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

ORFsvm_p2  <- ggplot(data = ORFsvm_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_ORFsvm <- ggarrange(ORFsvm_p1, ORFsvm_p2)
final_plot_ORFsvm <- annotate_figure(
  final_plot_ORFsvm,
  top = text_grob("ORFsvm", size = 15))
final_plot_ORFsvm

ROC_ORFsvm <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "ORFsvm",
                           metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_ORFsvm

predicciones_roc_ORFsvm <- predict(object = ROC_ORFsvm,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_ORFsvm <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_ORFsvm$TRUE.)

grafica_ROC_ORFsvm <- plot(curva_ROC_ORFsvm, main = "ORFsvm")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_ORFsvm, conf.level = 0.95)

predicciones_crudas_ORFsvm <- predict(ORFsvm_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_ORFsvm %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_ORFsvm, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_ORFsvm, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_ORFsvm != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Parallel Random Forest

parRF_df_am_final <- caret::train(Resultado ~ .,
                           data = df_am_entrenamiento,
                           method = "parRF",
                           trControl = fitControl_final,
                           tuneLenght = 40)

parRF_df_am_final$finalModel
parRF_df_am_final
parRF_df_am_final$resample %>% head(10)
summary(parRF_df_am_final$resample$Accuracy)

grafica_hiperparametros_parRF <- ggplot(parRF_df_am_final,  highlight = TRUE) +
  labs(title = "parRF") +
  theme_bw()

parRF_p1 <- ggplot(data = parRF_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(parRF_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

parRF_p2  <- ggplot(data = parRF_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_parRF <- ggarrange(parRF_p1, parRF_p2)
final_plot_parRF <- annotate_figure(
  final_plot_parRF,
  top = text_grob("parRF", size = 15))
final_plot_parRF


ROC_parRF <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "parRF",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_parRF

predicciones_roc_parRF <- predict(object = ROC_parRF,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_parRF <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_parRF$TRUE.)

grafica_ROC_parRF <- plot(curva_ROC_parRF, main = "parRF")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_parRF, conf.level = 0.95)

predicciones_crudas_parRF <- predict(parRF_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_parRF %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_parRF, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_parRF, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_parRF != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Random ferns

rFerns_df_am_final <- caret::train(Resultado ~ .,
                            data = df_am_entrenamiento,
                            method = "rFerns",
                            trControl = fitControl_final,
                            tuneLenght = 10)

rFerns_df_am_final$finalModel
rFerns_df_am_final
rFerns_df_am_final$resample %>% head(10)
summary(rFerns_df_am_final$resample$Accuracy)

grafica_hiperparametros_rFerns<- ggplot(rFerns_df_am_final,  highlight = TRUE) +
  labs(title = "rFerns") +
  theme_bw()

rFerns_p1 <- ggplot(data = rFerns_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(rFerns_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

rFerns_p2  <- ggplot(data = rFerns_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_rFerns <- ggarrange(rFerns_p1, rFerns_p2)
final_plot_rFerns <- annotate_figure(
  final_plot_rFerns,
  top = text_grob("rFerns", size = 15))
final_plot_rFerns

# No funciona ROC
# ROC_rFerns <- caret::train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "rFerns",
#                            metric = "ROC",
#                               trControl = ROC_control,
#                              family = "binomial")

# ROC_rFerns

# predicciones_roc_rFerns <- predict(object = ROC_rFerns,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_rFerns <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_rFerns$TRUE)

# plot(curva_ROC_rFerns)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_rFerns, conf.level = 0.95)

# predicciones_probabilidad_rFerns <- predict(rFerns_df_am_final, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_rFerns <- predict(rFerns_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_rFerns %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_rFerns, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_rFerns, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_rFerns != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Random Forest

ranger_df_am_final <- caret::train(Resultado ~ .,
                            data = df_am_entrenamiento,
                            method = "ranger",
                            trControl = fitControl_final,
                            tuneLenght = 40)

ranger_df_am_final$finalModel
ranger_df_am_final
ranger_df_am_final$resample %>% head(10)
summary(ranger_df_am_final$resample$Accuracy)

grafica_hiperparametros_ranger <- ggplot(ranger_df_am_final,  highlight = TRUE) +
  labs(title = "ranger") +
  theme_bw()

ranger_p1 <- ggplot(data = ranger_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(ranger_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

ranger_p2  <- ggplot(data = ranger_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_ranger <- ggarrange(ranger_p1, ranger_p2)
final_plot_ranger <- annotate_figure(
  final_plot_ranger,
  top = text_grob("ranger", size = 15))
final_plot_ranger

ROC_ranger <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "ranger",
                           metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

# ROC_ranger

predicciones_roc_ranger <- predict(object = ROC_ranger,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_ranger <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_ranger$TRUE.)

grafica_ROC_ranger <- plot(curva_ROC_ranger, main = "ranger")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_ranger, conf.level = 0.95)

predicciones_crudas_ranger <- predict(ranger_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_ranger %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_ranger, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_ranger, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_ranger != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")


# Random Forest

rf_df_am_final <- caret::train(Resultado ~ .,
                        data = df_am_entrenamiento,
                        method = "rf",
                        trControl = fitControl_final,
                        tuneLenght = 40)

rf_df_am_final$finalModel
rf_df_am_final
rf_df_am_final$resample %>% head(10)
summary(rf_df_am_final$resample$Accuracy)

grafica_hiperparametros_rf <- ggplot(rf_df_am_final,  highlight = TRUE) +
  labs(title = "precisión rf") +
  theme_bw()

rf_p1 <- ggplot(data = rf_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(rf_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

rf_p2  <- ggplot(data = rf_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_rf <- ggarrange(rf_p1, rf_p2)
final_plot_rf <- annotate_figure(
  final_plot_rf,
  top = text_grob("rf", size = 15))
final_plot_rf

ROC_rf <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "rf",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_rf

predicciones_roc_rf <- predict(object = ROC_rf,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_rf <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_rf$TRUE.)

grafica_ROC_rf <- plot(curva_ROC_rf, main = "rp")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_rf, conf.level = 0.95)

predicciones_crudas_rf <- predict(rf_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_rf %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_rf, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_rf, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_rf != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")
# Random Forest Rule-Based Model

rfRules_df_am_final <- caret::train(Resultado ~ .,
                             data = df_am_entrenamiento,
                             method = "rfRules",
                             trControl = fitControl_final,
                             tuneLenght = 40)

rfRules_df_am_final$finalModel
rfRules_df_am_final
rfRules_df_am_final$resample %>% head(10)
summary(rfRules_df_am_final$resample$Accuracy)

grafica_hiperparametros_rfRules <- ggplot(rfRules_df_am_final,  highlight = TRUE) +
  labs(title = "rfRules") +
  theme_bw()

rfRules_p1 <- ggplot(data = rfRules_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(rfRules_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

rfRules_p2  <- ggplot(data = rfRules_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_rfRules <- ggarrange(rfRules_p1, rfRules_p2)
final_plot_rfRules <- annotate_figure(
  final_plot_rfRules,
  top = text_grob("rfRules", size = 15))
final_plot_rfRules

# No funciona ROC
# ROC_rfRules <- caret::train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "rfRules",
#                               metric = "ROC",
#                               trControl = ROC_control,
#                              family = "binomial")

# ROC_rfRules

# predicciones_roc_rfRules <- predict(object = ROC_rfRules,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_rfRules <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_rfRules$TRUE)

# plot(curva_ROC_rfRules, main = "rfRules")

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_rfRules, conf.level = 0.95)

# predicciones_probabilidad_rfRules <- predict(rfRules_df_am_final, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_rfRules <- predict(rfRules_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_rfRules %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_rfRules, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_rfRules, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_rfRules != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Regularized Random Forest

RRF_df_am_final <- caret::train(Resultado ~ .,
                         data = df_am_entrenamiento,
                         method = "RRF",
                         trControl = fitControl_final,
                         tuneLenght = 40)

RRF_df_am_final$finalModel
RRF_df_am_final
RRF_df_am_final$resample %>% head(10)
summary(RRF_df_am_final$resample$Accuracy)

grafica_hiperparametros_RRF <- ggplot(RRF_df_am_final,  highlight = TRUE) +
  labs(title = "RRF") +
  theme_bw()

RRF_p1 <- ggplot(data = RRF_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(RRF_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

RRF_p2  <- ggplot(data = RRF_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_RRF <- ggarrange(RRF_p1, RRF_p2)
final_plot_RRF <- annotate_figure(
  final_plot_RRF,
  top = text_grob("RRF", size = 15))
final_plot_RRF


ROC_RRF <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "RRF",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_RRF

predicciones_roc_RRF <- predict(object = ROC_RRF,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_RRF <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_RRF$TRUE.)

grafica_ROC_RRF <- plot(curva_ROC_RRF, main = "RRF")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_RRF, conf.level = 0.95)

predicciones_crudas_RRF <- predict(RRF_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_RRF %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_RRF, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_RRF, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_RRF != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Regularized Random Forest

RRFglobal_df_am_final <- caret::train(Resultado ~ .,
                               data = df_am_entrenamiento,
                               method = "RRFglobal",
                               trControl = fitControl_final,
                               tuneLenght = 40)

RRFglobal_df_am_final$finalModel
RRFglobal_df_am_final
RRFglobal_df_am_final$resample %>% head(10)
summary(RRFglobal_df_am_final$resample$Accuracy)

grafica_hiperparametros_RRFglobal <- ggplot(RRFglobal_df_am_final,  highlight = TRUE) +
  labs(title = "RRFglobal") +
  theme_bw()

RRFglobal_p1 <- ggplot(data = RRFglobal_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(RRFglobal_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

RRFglobal_p2  <- ggplot(data = RRFglobal_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_RRFglobal <- ggarrange(RRFglobal_p1, RRFglobal_p2)
final_plot_RRFglobal <- annotate_figure(
  final_plot_RRFglobal,
  top = text_grob("RRFglobal", size = 15))
final_plot_RRFglobal


ROC_RRFglobal <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "RRFglobal",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_RRFglobal

predicciones_roc_RRFglobal <- predict(object = ROC_RRFglobal,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_RRFglobal <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_RRFglobal$TRUE.)

grafica_ROC_RRFglobal <- plot(curva_ROC_RRFglobal, main = "RRFglobal")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_RRFglobal, conf.level = 0.95)

predicciones_crudas_RRFglobal <- predict(RRFglobal_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_RRFglobal %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_RRFglobal, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_RRFglobal, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_RRFglobal != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Weighted Subspace Random Forest

wsrf_df_am_final <- caret::train(Resultado ~ .,
                          data = df_am_entrenamiento,
                          method = "wsrf",
                          trControl = fitControl_final,
                          tuneLenght = 40)

wsrf_df_am_final$finalModel
wsrf_df_am_final
wsrf_df_am_final$resample %>% head(10)
summary(wsrf_df_am_final$resample$Accuracy)

grafica_hiperparametros_wsrf <- ggplot(wsrf_df_am_final,  highlight = TRUE) +
  labs(title = "wsrf") +
  theme_bw()

wsrf_p1 <- ggplot(data = wsrf_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(wsrf_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

wsrf_p2  <- ggplot(data = wsrf_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_wsrf <- ggarrange(wsrf_p1, wsrf_p2)
final_plot_wsrf <- annotate_figure(
  final_plot_wsrf,
  top = text_grob("wsrf", size = 15))
final_plot_wsrf

ROC_wsrf <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "wsrf",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_wsrf

predicciones_roc_wsrf <- predict(object = ROC_wsrf,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_wsrf <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_wsrf$TRUE.)

grafica_ROC_wsrf <- plot(curva_ROC_wsrf, main = "wsrf")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_wsrf, conf.level = 0.95)

predicciones_crudas_wsrf <- predict(wsrf_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_wsrf %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_wsrf, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_wsrf, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_wsrf != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Rotation Forest

rotationForestCp_df_am_final <- caret::train(Resultado ~ .,
                                              data = df_am_entrenamiento,
                                              method = "rotationForestCp",
                                              trControl = fitControl_final,
                                              tuneLenght = 40)

rotationForestCp_df_am_final$finalModel
rotationForestCp_df_am_final
rotationForestCp_df_am_final$resample %>% head(10)
summary(rotationForestCp_df_am_final$resample$Accuracy)

grafica_hiperparametros_rotationForestCp <- ggplot(rotationForestCp_df_am_final,  highlight = TRUE) +
  labs(title = "rotationForestCp") +
  theme_bw()

rotationForestCp_p1 <- ggplot(data = rotationForestCp_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(rotationForestCp_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

rotationForestCp_p2  <- ggplot(data = rotationForestCp_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_rotationForestCp <- ggarrange(rotationForestCp_p1, rotationForestCp_p2)
final_plot_rotationForestCp <- annotate_figure(
  final_plot_rotationForestCp,
  top = text_grob("rotationForestCp", size = 15))
final_plot_rotationForestCp

ROC_rotationForestCp <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "rotationForestCp",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_rotationForestCp

predicciones_roc_rotationForestCp <- predict(object = ROC_rotationForestCp,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_rotationForestCp <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_rotationForestCp$TRUE.)

grafica_ROC_rotationForestCp <- plot(curva_ROC_rotationForestCp, main = "rotationForestCp")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_rotationForestCp, conf.level = 0.95)

predicciones_crudas_rotationForestCp <- predict(rotationForestCp_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_rotationForestCp %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_rotationForestCp, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_rotationForestCp, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_rotationForestCp != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Extreme Learning Machine")
# No esta disponible

elm_df_am_final <- caret::train(Resultado ~ .,
                                data = df_am_entrenamiento,
                                method = "elm", 
                                trControl = fitControl_final,
                                tuneLenght = 40)

elm_df_am_final$finalModel
elm_df_am_final
elm_df_am_final$resample %>% head(10)
summary(elm_df_am_final$resample$Accuracy)

grafica_hiperparametros_elm <- ggplot(elm_df_am_final,  highlight = TRUE) +
  labs(title = "elm") +
  theme_bw()

elm_p1 <- ggplot(data = elm_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(elm_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

elm_p2  <- ggplot(data = elm_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_elm <- ggarrange(elm_p1, elm_p2)
final_plot_elm <- annotate_figure(
  final_plot_elm,
  top = text_grob("elm", size = 15))
final_plot_elm

# No funciona ROC
# ROC_elm <- caret::train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "elm",
#                               metric = "ROC",
#                               trControl = ROC_control,
#                              family = "binomial")

# ROC_elm

# predicciones_roc_elm <- predict(object = ROC_elm,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_elm <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_elm$TRUE.)

# plot(curva_ROC_elm)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_elm, conf.level = 0.95)

# predicciones_probabilidad_elm <- predict(elm_df_am_final, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_elm <- predict(elm_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_elm %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_elm, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_elm, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_elm != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

## "Model Averaged Neural Network")

avNNet_df_am_final <- caret::train(Resultado ~ .,
                             data = df_am_entrenamiento,
                             method = "avNNet", 
                             trControl = fitControl_final,
                             tuneLenght = 40)

avNNet_df_am_final$finalModel
avNNet_df_am_final
avNNet_df_am_final$resample %>% head(10)
summary(avNNet_df_am_final$resample$Accuracy)

grafica_hiperparametros_avNNet <- ggplot(avNNet_df_am_final,  highlight = TRUE) +
  labs(title = "avNNet") +
  theme_bw()

avNNet_p1 <- ggplot(data = avNNet_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(avNNet_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

avNNet_p2  <- ggplot(data = avNNet_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_avNNet <- ggarrange(avNNet_p1, avNNet_p2)

final_plot_avNNet <- ggarrange(avNNet_p1, avNNet_p2)
final_plot_avNNet <- annotate_figure(
  final_plot_avNNet,
  top = text_grob("avNNet", size = 15))
final_plot_avNNet

ROC_avNNet <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "avNNet",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_avNNet

predicciones_roc_avNNet <- predict(object = ROC_avNNet,
                           newdata = df_am_test,
                           type = "prob")

curva_ROC_avNNet <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_avNNet$TRUE.)

grafica_ROC_avNNet <- plot(curva_ROC_avNNet, main = "avNNet")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_avNNet, conf.level = 0.95)

predicciones_crudas_avNNet <- predict(avNNet_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_avNNet %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_avNNet, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_avNNet, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_avNNet != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Neural Network

nnet_df_am_final <- caret::train(Resultado ~ .,
                           data = df_am_entrenamiento,
                           method = "nnet",
                           trControl = fitControl_final,
                           tuneLenght = 40)

nnet_df_am_final$finalModel
nnet_df_am_final
nnet_df_am_final$resample %>% head(10)
summary(nnet_df_am_final$resample$Accuracy)

grafica_hiperparametros_nnet <- ggplot(nnet_df_am_final,  highlight = TRUE) +
  labs(title = "nnet") +
  theme_bw()

nnet_p1 <- ggplot(data = nnet_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(nnet_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

nnet_p2  <- ggplot(data = nnet_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_nnet <- ggarrange(nnet_p1, nnet_p2)
final_plot_nnet <- annotate_figure(
  final_plot_nnet,
  top = text_grob("nnet", size = 15))
final_plot_nnet

ROC_nnet <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "nnet",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_nnet

predicciones_roc_nnet <- predict(object = ROC_nnet,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_nnet <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_nnet$TRUE.)

grafica_ROC_nnet <- plot(curva_ROC_nnet, main = "nnet")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_nnet, conf.level = 0.95)

predicciones_crudas_nnet <- predict(nnet_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_nnet %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_nnet, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_nnet, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_nnet != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Neural Networks with Feature Extraction

pcaNNet_df_am_final <- caret::train(Resultado ~ .,
                              data = df_am_entrenamiento,
                              method = "pcaNNet",
                              trControl = fitControl_final,
                              tuneLenght = 40)

pcaNNet_df_am_final$finalModel
pcaNNet_df_am_final
pcaNNet_df_am_final$resample %>% head(10)
summary(pcaNNet_df_am_final$resample$Accuracy)

grafica_hiperparametros_pcaNNet <- ggplot(pcaNNet_df_am_final,  highlight = TRUE) +
  labs(title = "pcaNNet") +
  theme_bw()

pcaNNet_p1 <- ggplot(data = pcaNNet_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(pcaNNet_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

pcaNNet_p2  <- ggplot(data = pcaNNet_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_pcaNNet <- ggarrange(pcaNNet_p1, pcaNNet_p2)
final_plot_pcaNNet <- annotate_figure(
  final_plot_pcaNNet,
  top = text_grob("pcaNNet", size = 15))
final_plot_pcaNNet

ROC_pcaNNet <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "pcaNNet",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_pcaNNet

predicciones_roc_pcaNNet <- predict(object = ROC_pcaNNet,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_pcaNNet <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_pcaNNet$TRUE.)

grafica_ROC_pcaNNet <- plot(curva_ROC_pcaNNet, main = "pcaNNet")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_pcaNNet, conf.level = 0.95)

predicciones_crudas_pcaNNet <- predict(pcaNNet_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_pcaNNet %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_pcaNNet, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_pcaNNet, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_pcaNNet != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Penalized Multinomial Regression

multinom_df_am_final <- caret::train(Resultado ~ .,
                               data = df_am_entrenamiento,
                               method = "multinom",
                               trControl = fitControl_final,
                               tuneLenght = 40)

multinom_df_am_final$finalModel
multinom_df_am_final
multinom_df_am_final$resample %>% head(10)
summary(multinom_df_am_final$resample$Accuracy)

grafica_hiperparametros_multinom <- ggplot(multinom_df_am_final,  highlight = TRUE) +
  labs(title = "multinom") +
  theme_bw()

multinom_p1 <- ggplot(data = multinom_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(multinom_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

multinom_p2  <- ggplot(data = multinom_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_multinom <- ggarrange(multinom_p1, multinom_p2)
final_plot_multinom <- annotate_figure(
  final_plot_multinom,
  top = text_grob("multinom", size = 15))
final_plot_multinom

ROC_multinom <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "multinom",
                              metric = "ROC",
                              trControl = ROC_control,
                              family = "binomial")

ROC_multinom

predicciones_roc_multinom <- predict(object = ROC_multinom,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_multinom <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_multinom$TRUE.)

grafica_ROC_multinom <- plot(curva_ROC_multinom, main = "multinom")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_multinom, conf.level = 0.95)

predicciones_crudas_multinom <- predict(multinom_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_multinom %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_multinom, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_multinom, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_multinom != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Radial Basis Function Network

rbfDDA_df_am_final <- caret::train(Resultado ~ .,
                             data = df_am_entrenamiento,
                             method = "rbfDDA",
                             trControl = fitControl_final,
                             tuneLenght = 40)

rbfDDA_df_am_final$finalModel
rbfDDA_df_am_final
rbfDDA_df_am_final$resample %>% head(10)
summary(rbfDDA_df_am_final$resample$Accuracy)

grafica_hiperparametros_rbfDDA <- ggplot(rbfDDA_df_am_final,  highlight = TRUE) +
  labs(title = "rbfDDA") +
  theme_bw()

rbfDDA_p1 <- ggplot(data = rbfDDA_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(rbfDDA_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

rbfDDA_p2  <- ggplot(data = rbfDDA_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_rbfDDA <- ggarrange(rbfDDA_p1, rbfDDA_p2)
final_plot_rbfDDA <- annotate_figure(
  final_plot_rbfDDA,
  top = text_grob("rbfDDA", size = 15))
final_plot_rbfDDA

ROC_rbfDDA <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "rbfDDA",
                              metric = "ROC",
                              trControl = ROC_control,
                              family = "binomial")

ROC_rbfDDA

predicciones_roc_rbfDDA <- predict(object = ROC_rbfDDA,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_rbfDDA <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_rbfDDA$TRUE.)

grafica_ROC_rbfDDA <- plot(curva_ROC_rbfDDA, main = "rbfDDA")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_rbfDDA, conf.level = 0.95)

predicciones_crudas_rbfDDA <- predict(rbfDDA_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_rbfDDA %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_rbfDDA, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_rbfDDA, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_rbfDDA != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Adjacent Categories Poobability Model for Ordinal Data

vglmAdjCat_df_am_final <- caret::train(Resultado ~ .,
                                 data = df_am_entrenamiento,
                                 method = "vglmAdjCat",
                                 trControl = fitControl_final,
                                 tuneLenght = 40)

vglmAdjCat_df_am_final$finalModel
vglmAdjCat_df_am_final
vglmAdjCat_df_am_final$resample %>% head(10)
summary(vglmAdjCat_df_am_final$resample$Accuracy)

grafica_hiperparametros_vglmAdjCat <- ggplot(vglmAdjCat_df_am_final,  highlight = TRUE) +
  labs(title = "vglmAdjCat") +
  theme_bw()

vglmAdjCat_p1 <- ggplot(data = vglmAdjCat_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(vglmAdjCat_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

vglmAdjCat_p2  <- ggplot(data = vglmAdjCat_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_vglmAdjCat <- ggarrange(vglmAdjCat_p1, vglmAdjCat_p2)
final_plot_vglmAdjCat <- annotate_figure(
  final_plot_vglmAdjCat,
  top = text_grob("vglmAdjCat", size = 15))
final_plot_vglmAdjCat

# No funciona ROC
# ROC_vglmAdjCat <- caret::train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "vglmAdjCat",
#                               metric = "ROC",
#                               trControl = ROC_control,
#                              family = "binomial")

# ROC_vglmAdjCat

# predicciones_roc_vglmAdjCat <- predict(object = ROC_vglmAdjCat,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_vglmAdjCat <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_vglmAdjCat$TRUE)

# plot(curva_ROC_vglmAdjCat)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_vglmAdjCat, conf.level = 0.95)

predicciones_crudas_vglmAdjCat <- predict(vglmAdjCat_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_vglmAdjCat %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_vglmAdjCat, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_vglmAdjCat, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_vglmAdjCat != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Boosted Logistic Regression

LogitBoost_df_am_final <- caret::train(Resultado ~ .,
                                 data = df_am_entrenamiento,
                                 method = "LogitBoost",
                                 trControl = fitControl_final,
                                 tuneLenght = 40)

LogitBoost_df_am_final$finalModel
LogitBoost_df_am_final
LogitBoost_df_am_final$resample %>% head(10)
summary(LogitBoost_df_am_final$resample$Accuracy)

grafica_hiperparametros_LogitBoost <- ggplot(LogitBoost_df_am_final,  highlight = TRUE) +
  labs(title = "LogitBoost") +
  theme_bw()

LogitBoost_p1 <- ggplot(data = LogitBoost_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(LogitBoost_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

LogitBoost_p2  <- ggplot(data = LogitBoost_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_LogitBoost <- ggarrange(LogitBoost_p1, LogitBoost_p2)
final_plot_LogitBoost <- annotate_figure(
  final_plot_LogitBoost,
  top = text_grob("LogitBoost", size = 15))
final_plot_LogitBoost

ROC_LogitBoost <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "LogitBoost",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_LogitBoost

predicciones_roc_LogitBoost <- predict(object = ROC_LogitBoost,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_LogitBoost <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_LogitBoost$TRUE.)

grafica_ROC_LogitBoost <- plot(curva_ROC_LogitBoost, main = "LogitBoost")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_LogitBoost, conf.level = 0.95)

# predicciones_probabilidad_LogitBoost <- predict(LogitBoost_df_am_final, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_LogitBoost <- predict(LogitBoost_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_LogitBoost %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_LogitBoost, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_LogitBoost, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_LogitBoost != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# # Continuation Ratio Model for Ordinal Data
# 
# vglmContRatio_df_am_final <- caret::train(Resultado ~ .,
#                                     data = df_am_entrenamiento,
#                                     method = "vglmContRatio",
#                                     trControl = fitControl_final,
#                                     tuneLenght = 5)
# 
# vglmContRatio_df_am_final$finalModel
# vglmContRatio_df_am_final
# vglmContRatio_df_am_final$resample %>% head(10)
# summary(vglmContRatio_df_am_final$resample$Accuracy)
# 
# grafica_hiperparametros_vglmContRatio <- ggplot(vglmContRatio_df_am_final,  highlight = TRUE) +
#   labs(title = "vglmContRatio") +
#   theme_bw()
# 
# vglmContRatio_p1 <- ggplot(data = vglmContRatio_df_am_final$resample, aes(x = Accuracy)) +
#   geom_density(alpha = 0.5, fill = "green") +
#   geom_vline(xintercept = mean(vglmContRatio_df_am_final$resample$Accuracy),
#              linetype = "dashed") +
#   labs(x = "precisión", y = "densidad")
# theme_bw()
# 
# vglmContRatio_p2  <- ggplot(data = vglmContRatio_df_am_final$resample, aes(x = 1, y = Accuracy)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
#   geom_jitter(width = 0.05) +
#   labs(y = "precisión") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# final_plot_vglmContRatio <- ggarrange(vglmContRatio_p1, vglmContRatio_p2)
# final_plot_vglmContRatio <- annotate_figure(
#   final_plot_vglmContRatio,
#   top = text_grob("vglmContRatio", size = 15))
# final_plot_vglmContRatio
# 
# # No funciona ROC
# # ROC_vglmContRatio <- caret::train(make.names(Resultado) ~ .,
# #                               data = df_am_entrenamiento,
# #                               method = "vglmContRatio",
# #                               metric = "ROC",
# #                               trControl = ROC_control,
# #                              family = "binomial")
# 
# # ROC_vglmContRatio
# 
# # predicciones_roc_vglmContRatio <- predict(object = ROC_vglmContRatio,
# #                            newdata = df_am_test,
# #                            type = "prob")
# # curva_ROC_vglmContRatio <- roc(response = df_am_test$Resultado,
# #                                      predictor = predicciones_roc_vglmContRatio$TRUE)
# 
# # plot(curva_ROC_vglmContRatio)
# 
# # Intervalo de confianza de la curva
# # ci.auc(curva_ROC_vglmContRatio, conf.level = 0.95)
# 
# # predicciones_probabilidad_vglmContRatio <- predict(vglmContRatio_df_am_final, newdata = df_am_test,
# #                                               type = "prob")
# 
# predicciones_crudas_vglmContRatio <- predict(vglmContRatio_df_am_final,
#                                                    newdata = df_am_test,
#                                                    type = "raw")
# predicciones_crudas_vglmContRatio %>% head(10)
# 
# caret::confusionMatrix(data = predicciones_crudas_vglmContRatio, reference = df_am_test$Resultado,
#                 positive = "TRUE")
# 
# CrossTable(predicciones_crudas_vglmContRatio, df_am_test$Resultado,
#            expected=FALSE, prop.r=TRUE, prop.c=TRUE,
#            prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
#            resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")
# 
# error_test <- mean(predicciones_crudas_vglmContRatio != df_am_test$Resultado)
# paste("Valor del error test: ", round(error_test*100, 2), "%")

# Cumulative Probability Model for Ordinal Data

vglmCumulative_df_am_final <- caret::train(Resultado ~ .,
                                           data = df_am_entrenamiento,
                                           method = "vglmCumulative",
                                           trControl = fitControl_final,
                                           tuneLenght = 40)

vglmCumulative_df_am_final$finalModel
vglmCumulative_df_am_final
vglmCumulative_df_am_final$resample %>% head(10)
summary(vglmCumulative_df_am_final$resample$Accuracy)

grafica_hiperparametros_vglmCumulative <- ggplot(vglmCumulative_df_am_final,  highlight = TRUE) +
  labs(title = "vglmCumulative") +
  theme_bw()

vglmCumulative_p1 <- ggplot(data = vglmCumulative_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(vglmCumulative_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

vglmCumulative_p2  <- ggplot(data = vglmCumulative_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_vglmCumulative <- ggarrange(vglmCumulative_p1, vglmCumulative_p2)
final_plot_vglmCumulative <- annotate_figure(
  final_plot_vglmCumulative,
  top = text_grob("vglmCumulative", size = 15))
final_plot_vglmCumulative

# No funciona ROC
# ROC_vglmCumulative <- caret::train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "vglmCumulative",
#                               metric = "ROC",
#                               trControl = ROC_control,
#                              family = "binomial")

# ROC_vglmCumulative

# predicciones_roc_vglmCumulative <- predict(object = ROC_vglmCumulative,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_vglmCumulative <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_vglmCumulative$TRUE)

# plot(curva_ROC_vglmCumulative)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_vglmCumulative, conf.level = 0.95)

# predicciones_probabilidad_vglmCumulative <- predict(vglmCumulative_df_am_final, newdata = df_am_test,
#                                               type = "prob")

predicciones_crudas_vglmCumulative <- predict(vglmCumulative_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_vglmCumulative %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_vglmCumulative, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_vglmCumulative, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_vglmCumulative != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Generalized Partial Least Squares

gpls_df_am_final <- caret::train(Resultado ~ .,
                           data = df_am_entrenamiento,
                           method = "gpls",
                           trControl = fitControl_final,
                           tuneLenght = 40)

gpls_df_am_final$finalModel
gpls_df_am_final
gpls_df_am_final$resample %>% head(10)
summary(gpls_df_am_final$resample$Accuracy)

grafica_hiperparametros_gpls <- ggplot(gpls_df_am_final,  highlight = TRUE) +
  labs(title = "gpls") +
  theme_bw()

gpls_p1 <- ggplot(data = gpls_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(gpls_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

gpls_p2  <- ggplot(data = gpls_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_gpls <- ggarrange(gpls_p1, gpls_p2)
final_plot_gpls <- annotate_figure(
  final_plot_gpls,
  top = text_grob("gpls", size = 15))
final_plot_gpls

ROC_gpls <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "gpls",
                              metric = "ROC",
                              trControl = ROC_control,
                              family = "binomial")

ROC_gpls

predicciones_roc_gpls <- predict(object = ROC_gpls,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_gpls <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_gpls$TRUE.)

grafica_ROC_gpls <- plot(curva_ROC_gpls, main = "gpls")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_gpls, conf.level = 0.95)

predicciones_crudas_gpls <- predict(gpls_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_gpls %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_gpls, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_gpls, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_gpls != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Regularized Logistic Regression

regLogistic_df_am_final <- caret::train(Resultado ~ .,
                                  data = df_am_entrenamiento,
                                  method = "regLogistic",
                                  trControl = fitControl_final,
                                  tuneLenght = 40)

regLogistic_df_am_final$finalModel
regLogistic_df_am_final
regLogistic_df_am_final$resample %>% head(10)
summary(regLogistic_df_am_final$resample$Accuracy)

grafica_hiperparametros_regLogistic <- ggplot(regLogistic_df_am_final,  highlight = TRUE) +
  labs(title = "regLogistic") +
  theme_bw()

regLogistic_p1 <- ggplot(data = regLogistic_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(regLogistic_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

regLogistic_p2  <- ggplot(data = regLogistic_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_regLogistic <- ggarrange(regLogistic_p1, regLogistic_p2)
final_plot_regLogistic <- annotate_figure(
  final_plot_regLogistic,
  top = text_grob("regLogistic", size = 15))
final_plot_regLogistic

ROC_regLogistic <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "regLogistic",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_regLogistic

predicciones_roc_regLogistic <- predict(object = ROC_regLogistic,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_regLogistic <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_regLogistic$TRUE.)

grafica_ROC_regLogistic <- plot(curva_ROC_regLogistic, main = "regLogistic")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_regLogistic, conf.level = 0.95)

predicciones_crudas_regLogistic <- predict(regLogistic_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_regLogistic %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_regLogistic, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_regLogistic, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_regLogistic != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# AdaBoost Clasification Trees

adaboost_df_am_final <- caret::train(Resultado ~ .,
                               data = df_am_entrenamiento,
                               method = "adaboost",
                               trControl = fitControl_final,
                               tuneLenght = 40)

adaboost_df_am_final$finalModel
adaboost_df_am_final
adaboost_df_am_final$resample %>% head(10)
summary(adaboost_df_am_final$resample$Accuracy)

grafica_hiperparametros_adaboost <- ggplot(adaboost_df_am_final,  highlight = TRUE) +
  labs(title = "adaboost") +
  theme_bw()

adaboost_p1 <- ggplot(data = adaboost_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(adaboost_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

adaboost_p2  <- ggplot(data = adaboost_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_adaboost <- ggarrange(adaboost_p1, adaboost_p2)
final_plot_adaboost <- annotate_figure(
  final_plot_adaboost,
  top = text_grob("adaboost", size = 15))
final_plot_adaboost

ROC_adaboost <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "adaboost",
                              metric = "ROC",
                              trControl = ROC_control,
                              family = "binomial")

ROC_adaboost

predicciones_roc_adaboost <- predict(object = ROC_adaboost,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_adaboost <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_adaboost$TRUE.)

grafica_ROC_adaboost <- plot(curva_ROC_adaboost, main = "adaboost")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_adaboost, conf.level = 0.95)

predicciones_crudas_adaboost <- predict(adaboost_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_adaboost %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_adaboost, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_adaboost, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_adaboost != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# AdaBoost.M1

AdaBoostM1_df_am_final <- caret::train(Resultado ~ .,
                                 data = df_am_entrenamiento,
                                 method = "AdaBoost.M1",
                                 trControl = fitControl_final,
                                 tuneLenght = 40)

AdaBoostM1_df_am_final$finalModel
AdaBoostM1_df_am_final
AdaBoostM1_df_am_final$resample %>% head(10)
summary(AdaBoostM1_df_am_final$resample$Accuracy)

grafica_hiperparametros_AdaBoostM1 <- ggplot(AdaBoostM1_df_am_final,  highlight = TRUE) +
  labs(title = "AdaBoostM1") +
  theme_bw()

AdaBoostM1_p1 <- ggplot(data = AdaBoostM1_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(AdaBoostM1_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

AdaBoostM1_p2  <- ggplot(data = AdaBoostM1_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_AdaBoostM1 <- ggarrange(AdaBoostM1_p1, AdaBoostM1_p2)
final_plot_AdaBoostM1 <- annotate_figure(
  final_plot_AdaBoostM1,
  top = text_grob("AdaBoostM1", size = 15))
final_plot_AdaBoostM1

ROC_AdaBoostM1 <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "AdaBoost.M1",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_AdaBoostM1

predicciones_roc_AdaBoostM1 <- predict(object = ROC_AdaBoostM1,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_AdaBoostM1 <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_AdaBoostM1$TRUE.)

grafica_ROC_AdaBoostM1 <- plot(curva_ROC_AdaBoostM1, main = "AdaBoostM1")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_AdaBoostM1, conf.level = 0.95)

predicciones_crudas_AdaBoostM1 <- predict(AdaBoostM1_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_AdaBoostM1 %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_AdaBoostM1, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_AdaBoostM1, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_AdaBoostM1 != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Bagged AdaBoost

AdaBag_df_am_final <- caret::train(Resultado ~ .,
                             data = df_am_entrenamiento,
                             method = "AdaBag",
                             trControl = fitControl_final,
                             tuneLenght = 40)

AdaBag_df_am_final$finalModel
AdaBag_df_am_final
AdaBag_df_am_final$resample %>% head(10)
summary(AdaBag_df_am_final$resample$Accuracy)

grafica_hiperparametros_AdaBag <- ggplot(AdaBag_df_am_final,  highlight = TRUE) +
  labs(title = "AdaBag") +
  theme_bw()

AdaBag_p1 <- ggplot(data = AdaBag_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(AdaBag_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

AdaBag_p2  <- ggplot(data = AdaBag_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_AdaBag <- ggarrange(AdaBag_p1, AdaBag_p2)
final_plot_AdaBag <- annotate_figure(
  final_plot_AdaBag,
  top = text_grob("AdaBag", size = 15))
final_plot_AdaBag

ROC_AdaBag <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "AdaBag",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_AdaBag

predicciones_roc_AdaBag <- predict(object = ROC_AdaBag,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_AdaBag <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_AdaBag$TRUE.)

grafica_ROC_AdaBag <- plot(curva_ROC_AdaBag, main = "AdaBag")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_AdaBag, conf.level = 0.95)

predicciones_crudas_AdaBag <- predict(AdaBag_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_AdaBag %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_AdaBag, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_AdaBag, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_AdaBag != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# K-nearest Neighbors

kknn_df_am_final <- caret::train(Resultado ~ .,
                           data = df_am_entrenamiento,
                           method = "kknn",
                           trControl = fitControl_final,
                           tuneLenght = 40)

kknn_df_am_final$finalModel
kknn_df_am_final
kknn_df_am_final$resample %>% head(10)
summary(kknn_df_am_final$resample$Accuracy)

grafica_hiperparametros_kknn <- ggplot(kknn_df_am_final,  highlight = TRUE) +
  labs(title = "kknn") +
  theme_bw()

kknn_p1 <- ggplot(data = kknn_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(kknn_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

kknn_p2  <- ggplot(data = kknn_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_kknn <- ggarrange(kknn_p1, kknn_p2)
final_plot_kknn <- annotate_figure(
  final_plot_kknn,
  top = text_grob("kknn", size = 15))
final_plot_kknn

ROC_kknn <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "kknn",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_kknn

predicciones_roc_kknn <- predict(object = ROC_kknn,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_kknn <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_kknn$TRUE.)

grafica_ROC_kknn <- plot(curva_ROC_kknn, main = "kknn")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_kknn, conf.level = 0.95)

predicciones_crudas_kknn <- predict(kknn_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_kknn %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_kknn, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_kknn, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_kknn != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# C5.0

C50_df_am_final <- caret::train(Resultado ~ .,
                          data = df_am_entrenamiento,
                          method = "C5.0",
                          trControl = fitControl_final,
                          tuneLenght = 40)

# C50_df_am_final$finalModel
C50_df_am_final
C50_df_am_final$resample %>% head(10)
summary(C50_df_am_final$resample$Accuracy)

grafica_hiperparametros_C50 <- ggplot(C50_df_am_final,  highlight = TRUE) +
  labs(title = "C50") +
  theme_bw()

C50_p1 <- ggplot(data = C50_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(C50_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

C50_p2  <- ggplot(data = C50_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_C50 <- ggarrange(C50_p1, C50_p2)
final_plot_C50 <- annotate_figure(
  final_plot_C50,
  top = text_grob("C50", size = 15))
final_plot_C50

ROC_C50 <- caret::train(make.names(Resultado) ~ .,
                              data = df_am_entrenamiento,
                              method = "C5.0",
                              metric = "ROC",
                              trControl = ROC_control,
                             family = "binomial")

ROC_C50

predicciones_roc_C50 <- predict(object = ROC_C50,
                           newdata = df_am_test,
                           type = "prob")
curva_ROC_C50 <- roc(response = df_am_test$Resultado,
                                     predictor = predicciones_roc_C50$TRUE.)

grafica_ROC_C50 <- plot(curva_ROC_C50, main = "C50")

# Intervalo de confianza de la curva
ci.auc(curva_ROC_C50, conf.level = 0.95)

predicciones_crudas_C50 <- predict(C50_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_C50 %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_C50, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_C50, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_C50 != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Cost-Sensitive C5.0

C50Cost_df_am_final <- caret::train(Resultado ~ .,
                              data = df_am_entrenamiento,
                              method = "C5.0Cost",
                              trControl = fitControl_final,
                              tuneLenght = 40)

#C50Cost_df_am_final$finalModel
C50Cost_df_am_final
C50Cost_df_am_final$resample %>% head(10)
summary(C50Cost_df_am_final$resample$Accuracy)

grafica_hiperparametros_C50Cost <- ggplot(C50Cost_df_am_final,  highlight = TRUE) +
  labs(title = "C50Cost") +
  theme_bw()

C50Cost_p1 <- ggplot(data = C50Cost_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(C50Cost_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

C50Cost_p2  <- ggplot(data = C50Cost_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_C50Cost <- ggarrange(C50Cost_p1, C50Cost_p2)
final_plot_C50Cost <- annotate_figure(
  final_plot_C50Cost,
  top = text_grob("C50Cost", size = 15))
final_plot_C50Cost

# No funciona ROC
# ROC_C50Cost <- caret::train(make.names(Resultado) ~ .,
#                         data = df_am_entrenamiento,
#                         method = "C5.0Cost",
#                         metric = "ROC",
#                         trControl = ROC_control,
#                         family = "binomial")
# 
# ROC_C50Cost
# 
# predicciones_roc_C50Cost <- predict(object = ROC_C50Cost,
#                                 newdata = df_am_test,
#                                 type = "prob")
# curva_ROC_C50Cost <- roc(response = df_am_test$Resultado,
#                      predictor = predicciones_roc_C50Cost$TRUE.)
# 
# plot(curva_ROC_C50Cost)
# 
# # Intervalo de confianza de la curva
# ci.auc(curva_ROC_C50, conf.level = 0.95)

predicciones_crudas_C50Cost <- predict(C50Cost_df_am_final,
                                   newdata = df_am_test,
                                   type = "raw")

predicciones_crudas_C50Cost %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_C50Cost, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_C50Cost, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_C50Cost != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Localized Linear Discriminant Analysis

# loclda_df_am_final <- caret::train(Resultado ~ .,
#                              data = df_am_entrenamiento,
#                              method = "loclda",
#                              trControl = fitControl_final,
#                              tuneLenght = 5)
# 
# loclda_df_am_final$finalModel
# loclda_df_am_final
# loclda_df_am_final$resample %>% head(10)
# summary(loclda_df_am_final$resample$Accuracy)
# 
# grafica_hiperparametros_loclda <- ggplot(loclda_df_am_final,  highlight = TRUE) +
#   labs(title = "loclda") +
#   theme_bw()
# 
# loclda_p1 <- ggplot(data = loclda_df_am_final$resample, aes(x = Accuracy)) +
#   geom_density(alpha = 0.5, fill = "green") +
#   geom_vline(xintercept = mean(loclda_df_am_final$resample$Accuracy),
#              linetype = "dashed") +
#   labs(x = "precisión", y = "densidad")
# theme_bw()
# 
# loclda_p2  <- ggplot(data = loclda_df_am_final$resample, aes(x = 1, y = Accuracy)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
#   geom_jitter(width = 0.05) +
#   labs(y = "precisión") +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# final_plot_loclda <- ggarrange(loclda_p1, loclda_p2)
# final_plot_loclda <- annotate_figure(
#   final_plot_loclda,
#   top = text_grob("loclda", size = 15))
# final_plot_loclda
# 
# # No funciona ROC
# # ROC_loclda <- caret::train(make.names(Resultado) ~ .,
# #                               data = df_am_entrenamiento,
# #                               method = "loclda",
# #                               metric = "ROC",
# #                               trControl = ROC_control,
# #                              family = "binomial")
# 
# # ROC_loclda
# 
# # predicciones_roc_loclda <- predict(object = ROC_loclda,
# #                            newdata = df_am_test,
# #                            type = "prob")
# # curva_ROC_loclda <- roc(response = df_am_test$Resultado,
# #                                      predictor = predicciones_roc_loclda$TRUE)
# 
# # plot(curva_ROC_loclda)
# 
# # Intervalo de confianza de la curva
# # ci.auc(curva_ROC_loclda, conf.level = 0.95)
# 
# # predicciones_probabilidad_loclda <- predict(loclda_df_am_final, newdata = df_am_test,
# #                                               type = "prob")
# 
# predicciones_crudas_loclda <- predict(loclda_df_am_final,
#                                                    newdata = df_am_test,
#                                                    type = "raw")
# predicciones_crudas_loclda %>% head(10)
# 
# caret::confusionMatrix(data = predicciones_crudas_loclda, reference = df_am_test$Resultado,
#                 positive = "TRUE")
# 
# CrossTable(predicciones_crudas_loclda, df_am_test$Resultado,
#            expected=FALSE, prop.r=TRUE, prop.c=TRUE,
#            prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
#            resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")
# 
# error_test <- mean(predicciones_crudas_loclda != df_am_test$Resultado)
# paste("Valor del error test: ", round(error_test*100, 2), "%")

# Optimal Weighted Nearest Neighbor Classifier

ownn_df_am_final <- caret::train(Resultado ~ .,
                           data = df_am_entrenamiento,
                           method = "ownn",
                           trControl = fitControl_final,
                           tuneLenght = 40)

ownn_df_am_final$finalModel
ownn_df_am_final
ownn_df_am_final$resample %>% head(10)
summary(ownn_df_am_final$resample$Accuracy)

grafica_hiperparametros_ownn <- ggplot(ownn_df_am_final,  highlight = TRUE) +
  labs(title = "ownn") +
  theme_bw()

ownn_p1 <- ggplot(data = ownn_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(ownn_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

ownn_p2  <- ggplot(data = ownn_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_ownn <- ggarrange(ownn_p1, ownn_p2)
final_plot_ownn <- annotate_figure(
  final_plot_ownn,
  top = text_grob("ownn", size = 15))
final_plot_ownn

# No funciona ROC
# ROC_ownn <- caret::train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "ownn",
#                               metric = "ROC",
#                               trControl = ROC_control,
#                              family = "binomial")

# ROC_ownn

# predicciones_roc_ownn <- predict(object = ROC_ownn,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_ownn <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_ownn$TRUE)

# plot(curva_ROC_ownn)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_ownn, conf.level = 0.95)

predicciones_crudas_ownn <- predict(ownn_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_ownn %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_ownn, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_ownn, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_ownn != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

# Stabilized nearest Neighbor Classifier

snn_df_am_final <- caret::train(Resultado ~ .,
                          data = df_am_entrenamiento,
                          method = "snn",
                          trControl = fitControl_final,
                          tuneLenght = 40)

snn_df_am_final$finalModel
snn_df_am_final
snn_df_am_final$resample %>% head(10)
summary(snn_df_am_final$resample$Accuracy)

grafica_hiperparametros_snn <- ggplot(snn_df_am_final,  highlight = TRUE) +
  labs(title = "snn") +
  theme_bw()

snn_p1 <- ggplot(data = snn_df_am_final$resample, aes(x = Accuracy)) +
  geom_density(alpha = 0.5, fill = "green") +
  geom_vline(xintercept = mean(snn_df_am_final$resample$Accuracy),
             linetype = "dashed") +
  labs(x = "precisión", y = "densidad")
theme_bw()

snn_p2  <- ggplot(data = snn_df_am_final$resample, aes(x = 1, y = Accuracy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "green") +
  geom_jitter(width = 0.05) +
  labs(y = "precisión") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

final_plot_snn <- ggarrange(snn_p1, snn_p2)
final_plot_snn <- annotate_figure(
  final_plot_snn,
  top = text_grob("snn", size = 15))
final_plot_snn

# No funciona ROC
# ROC_snn <- caret::train(make.names(Resultado) ~ .,
#                               data = df_am_entrenamiento,
#                               method = "snn",
#                               metric = "ROC",
#                               trControl = ROC_control,
#                              family = "binomial")

# ROC_snn

# predicciones_roc_snn <- predict(object = ROC_snn,
#                            newdata = df_am_test,
#                            type = "prob")
# curva_ROC_snn <- roc(response = df_am_test$Resultado,
#                                      predictor = predicciones_roc_snn$TRUE)

# plot(curva_ROC_snn)

# Intervalo de confianza de la curva
# ci.auc(curva_ROC_snn, conf.level = 0.95)

predicciones_crudas_snn <- predict(snn_df_am_final,
                                                   newdata = df_am_test,
                                                   type = "raw")
predicciones_crudas_snn %>% head(10)

caret::confusionMatrix(data = predicciones_crudas_snn, reference = df_am_test$Resultado,
                positive = "TRUE")

CrossTable(predicciones_crudas_snn, df_am_test$Resultado,
           expected=FALSE, prop.r=TRUE, prop.c=TRUE,
           prop.t=TRUE,prop.chisq=TRUE, chisq = TRUE, fisher=TRUE,  
           resid=TRUE, sresid=TRUE, asresid=TRUE, format = "SPSS")

error_test <- mean(predicciones_crudas_snn != df_am_test$Resultado)
paste("Valor del error test: ", round(error_test*100, 2), "%")

#### Comparacion de modelos

modelos_am <- list(svmLinearWeights2 = svmLinearWeights2_df_am_final,
                svmLinearWeights = svmLinearWeights_df_am_final,
                svmRadialWeights = svmRadialWeights_df_am_final,
                svmLinear3 = svmLinear3_df_am_final,
                # lssvmPoly = lssvmPoly_final,
                lssvmRadial = lssvmRadial_df_am_final,
                #svmlinear = svmLinear_df_am_final,
                svmLinear2 = svmLinear2_df_am_final,
                svmPoly = svmPoly_df_am_final,
                svmRadial = svmRadial_df_am_final,
                svmRadialCost = svmRadialCost_df_am_final,
                svmRadialSigma = svmRadialSigma_df_am_final,
                mlp = mlp_df_am_final,
                mlpWeightDecay = mlpWeightDecay_df_am_final,
                mlpWeightDecayML = mlpWeightDecayML_df_am_final,
                mplML = mlpML_df_am_final,
                ORFlog = ORFlog_df_am_final,
                ORFpls = ORFpls_df_am_final,
                ORFridge = ORFridge_df_am_final,
                ORFsvm = ORFsvm_df_am_final,
                parRF = parRF_df_am_final,
                rFerns = rFerns_df_am_final,
                ranger = ranger_df_am_final,
                rf = rf_df_am_final,
                rfRules = rfRules_df_am_final,
                RRF = RRF_df_am_final,
                RRFglobal = RRFglobal_df_am_final,
                wsrf = wsrf_df_am_final,
                rotationForestCp = rotationForestCp_df_am_final,
                elm = elm_df_am_final,
                avNNet = avNNet_df_am_final,
                nnet = nnet_df_am_final,
                pcaNNet = pcaNNet_df_am_final,
                multinom = multinom_df_am_final,
                rbfDDA = rbfDDA_df_am_final,
                vglmAdjCat = vglmAdjCat_df_am_final,
                LogitBoost = LogitBoost_df_am_final,
                # vglmContRatio = vglmContRatio_df_am_final,
                vglmCumulative = vglmCumulative_df_am_final,
                gpls = gpls_df_am_final,
                multinom = multinom_df_am_final,
                regLogistic = regLogistic_df_am_final,
                adaboost = adaboost_df_am_final,
                Adaboost = AdaBoostM1_df_am_final,
                AdaBag = AdaBag_df_am_final,
                kknn = kknn_df_am_final,
                C50 = C50_df_am_final,
                C50Cost = C50Cost_df_am_final,
                # loclda = loclda_df_am_final,
                ownn = ownn_df_am_final,
                snn = snn_df_am_final
                # ordinalRF = ordinalRF_df_am,
                # monmlp = monmlp_df_am_final,
                # cforest = cforest_df_am,
                # mlpKerasDropoutCost = mlpKerasDropoutCost_df_am,
                # mlpKerasDropout = mlpKerasDropout_df_am,
                # mlpKerasDecayCost = mlpKerasDecayCost_df_am,
                # Rborist = Rborist_df_am,
                # rotationForest = rotationForest_df_am,
                # dnn = dnn_df_am,
                # bayesglm = bayesglm_df_am,
                # LMT = LMT_df_am,
                # plr = plr_df_am,
                # ada = ada_df_am,
                # glmboost = glmboost_df_am,
                # BstLm = BstLm_df_am,
                # blackboost = blackboost_df_am,
                # bstTree = bstTree_df_am,
                # deepboost = deepboost_df_am,
                # xgbDART = xgbDART_df_am,
                # xgbLinear = xgbLinear_df_am,
                # gbm = gbm_df_am,
                # knn = knn_df_am,
                # extraTrees = extraTrees_df_am,
                # xgbTree = xgbTree_df_am
                )

resultados_modelos_resamples_am <- resamples(modelos_am)
resultados_modelos_resamples_am$values %>% head(10)

summary(resultados_modelos_resamples_am)

metricas_modelos_resamples_am <- resultados_modelos_resamples_am$values %>%
                         gather(key = "modelo", value = "valor", -Resample) %>%
                         separate(col = "modelo", into = c("modelo", "metrica"),
                                  sep = "~", remove = TRUE)
metricas_modelos_resamples_am %>% head()


metricas_modelos_resamples_am %>%
  group_by(modelo, metrica) %>% 
  summarise(media = mean(valor)) %>%
  spread(key = metrica, value = media) %>%
  arrange(desc(Accuracy)) %>% print(n=50)

metricas_modelos_resamples_am %>%
  filter(metrica == "Accuracy") %>%
  group_by(modelo) %>%
  summarise(media = mean(valor)) %>%
  ggplot(aes(x = reorder(modelo, media), y = media, label = round(media, 2))) +
    geom_segment(aes(x = reorder(modelo, media), y = 0,
                     xend = modelo, yend = media),
                     color = "grey50") +
    geom_point(size = 7, color = "green") +
    geom_text(color = "black", size = 2.65) +
    scale_y_continuous(limits = c(0, 1)) +
    # Precisión a nivel basal
    geom_hline(yintercept = 0.72, linetype = "dashed") +
    annotate(geom = "text", y = 0.72, x = 10, label = "Precisión basal") +
    labs(title = "Precisión promedio por adaptive-CV",
         subtitle = "Modelos ordenados por media",
         x = "modelo") +
    coord_flip() +
    theme_bw()

metricas_modelos_resamples_am %>%
  filter(metrica == "Accuracy") %>%
  group_by(modelo) %>%
  mutate(media = mean(valor)) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(modelo, media), y = valor, color = modelo)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.6) +
    scale_y_continuous(limits = c(0, 1)) +
    # Precisión a nivel basal
    geom_hline(yintercept = 0.72, linetype = "dashed") +
    annotate(geom = "text", y = 0.72, x = 8.5, label = "Precisión basal") +
    theme_bw() +
    labs(title = "Precisión media por adaptative-CV",
         subtitle = "Modelos ordenados por media") +
    coord_flip() +
    theme(legend.position = "none")

# Test de wilcoxon

metricas_precisión_am <- metricas_modelos_resamples_am %>% filter(metrica == "Accuracy")
comparaciones_am  <- pairwise.wilcox.test(x = metricas_precisión_am$valor,
                                        g = metricas_precisión_am$modelo,
                                        paired = TRUE,
                                        p.adjust.method = "holm")

# comparaciones_am <- comparaciones_am$p.value %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "modeloA") %>%
#   gather(key = "modeloB", value = "p_value", -modeloA) %>%
#   na.omit() %>%
#   arrange(modeloA)
# 
# comparaciones_am

# Error de Test
predicciones_am <- extractPrediction(
                  models = modelos_am,
                  testX = df_am_test[, !names(df_am_test) %in% c("Resultado")],
                  testY = df_am_test$Resultado
                )
predicciones_am %>% head()

metricas_predicciones_am <- predicciones_am %>%
                         mutate(acierto = ifelse(obs == pred, TRUE, FALSE)) %>%
                         group_by(object, dataType) %>%
                         summarise(accuracy = mean(acierto))

metricas_predicciones_am %>%
  spread(key = dataType, value = accuracy) %>%
  arrange(desc(Test))

ggplot(data = metricas_predicciones_am,
       aes(x = reorder(object, accuracy), y = accuracy,
           color = dataType, label = round(accuracy, 2))) +
  geom_point(size = 7) +
  scale_color_manual(values = c("red", "green")) +
  geom_text(color = "black", size = 2.65) +
  scale_y_continuous(limits = c(0, 1)) +
  # Accuracy basal
  geom_hline(yintercept = 0.72, linetype = "dashed") +
  annotate(geom = "text", y = 0.72, x = 8.5, label = "Precisión basal") +
  coord_flip() +
  labs(title = "Precisión de entrenamiento y test",
       x = "modelo",
       y = "precisión") +
  theme_bw() +
  theme(legend.position = "bottom")

#Graficas de hiperametros unificadas

lista_grafica_final_hiperparametros <-list(grafica_hiperparametros_svmLinearWeights2, 
                               grafica_hiperparametros_svmLinearWeights,
                               grafica_hiperparametros_svmRadialWeights,
                               grafica_hiperparametros_svmLinear3,
                               #grafica_hiperparametros_lssvmPoly,
                               grafica_hiperparametros_lssvmRadial,
                               # grafica_hiperparametros_svmLinear,
                               grafica_hiperparametros_svmLinear2,
                               grafica_hiperparametros_svmPoly,
                               grafica_hiperparametros_svmRadial,
                               grafica_hiperparametros_svmRadialCost,
                               grafica_hiperparametros_svmRadialSigma,
                               grafica_hiperparametros_mlp,
                               grafica_hiperparametros_mlpWeightDecay,
                               grafica_hiperparametros_mlpWeightDecayML,
                               grafica_hiperparametros_mlpML,
                               grafica_hiperparametros_ORFlog,
                               grafica_hiperparametros_ORFpls,
                               grafica_hiperparametros_ORFridge,
                               grafica_hiperparametros_ORFsvm,
                               grafica_hiperparametros_parRF,
                               grafica_hiperparametros_rFerns,
                               grafica_hiperparametros_ranger,
                               grafica_hiperparametros_rf,
                               grafica_hiperparametros_rfRules,
                               grafica_hiperparametros_RRF,
                               grafica_hiperparametros_RRFglobal,
                               grafica_hiperparametros_wsrf,
                               grafica_hiperparametros_rotationForestCp,
                               grafica_hiperparametros_elm,
                               grafica_hiperparametros_avNNet,
                               grafica_hiperparametros_nnet,
                               grafica_hiperparametros_pcaNNet,
                               grafica_hiperparametros_multinom,
                               grafica_hiperparametros_rbfDDA,
                               grafica_hiperparametros_vglmAdjCat,
                               grafica_hiperparametros_LogitBoost,
                               # grafica_hiperparametros_vglmContRatio,
                               grafica_hiperparametros_vglmCumulative,
                               grafica_hiperparametros_gpls,
                               grafica_hiperparametros_multinom,
                               grafica_hiperparametros_regLogistic,
                               grafica_hiperparametros_adaboost,
                               grafica_hiperparametros_AdaBoostM1,
                               grafica_hiperparametros_AdaBag,
                               grafica_hiperparametros_kknn,
                               grafica_hiperparametros_C50,
                               grafica_hiperparametros_C50Cost,
                               #grafica_hiperparametros_loclda,
                               grafica_hiperparametros_ownn,
                               grafica_hiperparametros_snn) 

Grafica_final_hiperparametros1 <- ggarrange(plotlist = lista_grafica_final_hiperparametros[1:4])
Grafica_final_hiperparametros1 <- annotate_figure(Grafica_final_hiperparametros1, top = text_grob("Precisión de hiperparametros", size = 20))

Grafica_final_hiperparametros2 <- ggarrange(plotlist = lista_grafica_final_hiperparametros[5:8])

Grafica_final_hiperparametros3 <- ggarrange(plotlist = lista_grafica_final_hiperparametros[9:17], ncol = 3, nrow = 3)

Grafica_final_hiperparametros4 <- ggarrange(plotlist = lista_grafica_final_hiperparametros[18:26], ncol = 3, nrow = 3)

Grafica_final_hiperparametros5 <- ggarrange(plotlist = lista_grafica_final_hiperparametros[27:35], ncol = 3, nrow = 3)

Grafica_final_hiperparametros6 <- ggarrange(plotlist = lista_grafica_final_hiperparametros[36:44], ncol = 3, nrow = 3)

Grafica_final_hiperparametros7 <- ggarrange(plotlist = lista_grafica_final_hiperparametros[45:49])

Grafica_final_hiperparametros1
Grafica_final_hiperparametros2
Grafica_final_hiperparametros3
Grafica_final_hiperparametros4
Grafica_final_hiperparametros5
Grafica_final_hiperparametros6 
Grafica_final_hiperparametros7 
#Graficas de precisión unificadas

lista_grafica_final_precision <-list(final_plot_svmLinearWeights2, 
                                     final_plot_svmLinearWeights,
                                     final_plot_svmRadialWeights,
                                     final_plot_svmLinear3,
                                     #final_plot_lssvmPoly,
                                     final_plot_lssvmRadial,
                                     # final_plot_svmLinear,
                                     final_plot_svmLinear2,
                                     final_plot_svmPoly,
                                     final_plot_svmRadial,
                                     final_plot_svmRadialCost,
                                     final_plot_svmRadialSigma,
                                     final_plot_mlp,
                                     final_plot_mlpWeightDecay,
                                     final_plot_mlpWeightDecayML,
                                     final_plot_mlpML,
                                     final_plot_ORFlog,
                                     final_plot_ORFpls,
                                     final_plot_ORFridge,
                                     final_plot_ORFsvm,
                                     final_plot_parRF,
                                     final_plot_rFerns,
                                     final_plot_ranger,
                                     final_plot_rf,
                                     final_plot_rfRules,
                                     final_plot_RRF,
                                     final_plot_RRFglobal,
                                     final_plot_wsrf,
                                     final_plot_rotationForestCp,
                                     final_plot_elm,
                                     final_plot_avNNet,
                                     final_plot_nnet,
                                     final_plot_pcaNNet,
                                     final_plot_multinom,
                                     final_plot_rbfDDA,
                                     final_plot_vglmAdjCat,
                                     final_plot_LogitBoost,
                                     # final_plot_vglmContRatio,
                                     final_plot_vglmCumulative,
                                     final_plot_gpls,
                                     final_plot_multinom,
                                     final_plot_regLogistic,
                                     final_plot_adaboost,
                                     final_plot_AdaBoostM1,
                                     final_plot_AdaBag,
                                     final_plot_kknn,
                                     final_plot_C50,
                                     final_plot_C50Cost,
                                     # final_plot_loclda,
                                     final_plot_ownn,
                                     final_plot_snn) 

Grafica_final_precision1 <- ggarrange(plotlist = lista_grafica_final_precision[1:12], ncol = 4, nrow = 3)
Grafica_final_precision1 <- annotate_figure(Grafica_final_precision1, top = text_grob("Precisión de cada modelo", size = 20))

Grafica_final_precision2 <- ggarrange(plotlist = lista_grafica_final_precision[13:24], ncol = 4, nrow = 3)

Grafica_final_precision3 <- ggarrange(plotlist = lista_grafica_final_precision[25:36], ncol = 4, nrow = 3)

Grafica_final_precision4 <- ggarrange(plotlist = lista_grafica_final_precision[37:45], ncol = 3, nrow = 3)

Grafica_final_precision5 <- ggarrange(plotlist = lista_grafica_final_precision[46:49])

Grafica_final_precision1
Grafica_final_precision2
Grafica_final_precision3
Grafica_final_precision4
Grafica_final_precision5

# Graficas ROC unificadas
# ggarrange no funciona para objetos roc
# 
# lista_grafica_final_ROC <-list(# grafica_ROC_svmLinearWeights2, 
#                                      grafica_ROC_svmLinearWeights,
#                                      grafica_ROC_svmRadialWeights,
#                                      # grafica_ROC_svmLinear3,
#                                      # grafica_ROC_lssvmPoly,
#                                      # grafica_ROC_lssvmRadial,
#                                      grafica_ROC_svmLinear,
#                                      grafica_ROC_svmLinear2,
#                                      grafica_ROC_svmPoly,
#                                      grafica_ROC_svmRadial,
#                                      grafica_ROC_svmRadialCost,
#                                      grafica_ROC_svmRadialSigma,
#                                      grafica_ROC_mlp,
#                                      grafica_ROC_mlpWeightDecay,
#                                      grafica_ROC_mlpWeightDecayML,
#                                      grafica_ROC_mlpML,
#                                      grafica_ROC_ORFlog,
#                                      grafica_ROC_ORFpls,
#                                      grafica_ROC_ORFridge,
#                                      grafica_ROC_ORFsvm,
#                                      grafica_ROC_parRF,
#                                      # grafica_ROC_rFerns,
#                                      grafica_ROC_ranger,
#                                      grafica_ROC_rf,
#                                      # grafica_ROC_rfRules,
#                                      grafica_ROC_RRF,
#                                      grafica_ROC_RRFglobal,
#                                      grafica_ROC_wsrf,
#                                      grafica_ROC_rotationForestCp,
#                                      # grafica_ROC_elm,
#                                      grafica_ROC_avNNet,
#                                      grafica_ROC_nnet,
#                                      grafica_ROC_pcaNNet,
#                                      grafica_ROC_multinom,
#                                      grafica_ROC_rbfDDA,
#                                      # grafica_ROC_vglmAdjCat,
#                                      grafica_ROC_LogitBoost,
#                                      # grafica_ROC_vglmContRatio,
#                                      # grafica_ROC_vglmCumulative,
#                                      grafica_ROC_gpls,
#                                      grafica_ROC_multinom,
#                                      grafica_ROC_regLogistic,
#                                      grafica_ROC_adaboost,
#                                      grafica_ROC_AdaBoostM1,
#                                      grafica_ROC_AdaBag,
#                                      grafica_ROC_kknn,
#                                      grafica_ROC_C50
#                                      # grafica_ROC_C50Cost,
#                                      # grafica_ROC_loclda,
#                                      # grafica_ROC_ownn,
#                                      # grafica_ROC_snn
#                                      )
# 
# Grafica_final_ROC1 <- ggarrange(plotlist = lista_grafica_final_ROC[1:12], ncol = 4, nrow = 3)
# 
# Grafica_final_ROC1 <- annotate_figure(Grafica_final_ROC1, top = text_grob("Graficas ROC", size = 20))
# 
# Grafica_final_ROC2 <- ggarrange(plotlist = lista_grafica_final_ROC[13:24], ncol = 4, nrow = 3)
# Grafica_final_ROC2 <- annotate_figure(Grafica_final_ROC2, top = text_grob("Graficas ROC", size = 20))
# 
# Grafica_final_ROC3 <- ggarrange(plotlist = lista_grafica_final_ROC[25:36], ncol = 4, nrow = 3)
# Grafica_final_ROC3 <- annotate_figure(Grafica_final_ROC3, top = text_grob("Graficas ROC", size = 20))
# 
# Grafica_final_ROC4 <- ggarrange(plotlist = lista_grafica_final_ROC[37:37])
# Grafica_final_ROC4 <- annotate_figure(Grafica_final_ROC4, top = text_grob("Graficas ROC", size = 20))
# 
# Grafica_final_ROC1
# Grafica_final_ROC2
# Grafica_final_ROC3
# Grafica_final_ROC4
# 
# 
sink(type="output")
sink(type="message")
close(output)  # Close connection to log file

# comando para guardar en el disco el modelo elegido
# saveRDS(nombre modelo, 'modelo.rds')

saveRDS(snn_df_am_final, "modelo_am.rds")
saveRDS(trained_recipe, "trained_recipe_am.rds")
saveRDS(rf_ga_am$optVariables, "predictores_filtrados_am.rds")

## Validación final
# y_preds_cart <- predict(cart_p21, df_test)
# sqrt(mean((y_preds_cart - df_test$p21)**2))