##########
# Charging all the csv created for each dataset. 

GSE97362 <- read.csv( "path/to/the/sample/sheet/created/for/GSE/97362_sample_sheet.csv")

############
# Load the dplyr package

{library(dplyr)
library(purrr)
library(stringr)}

# Function to separate the data into disease states

Sep <- split(GSE97362, GSE97362$sample.type.ch1)

## Grouping by type of samples

Kabuki <- rbind(Sep$`KMT2D LOF discovery cohort`,Sep$`Control for KMT2D LOF discovery cohort`)
Control <- Sep$`Control for validation cohort`

rm(Sep, GSE97362)

# Making sure that age and Sex are present 

Control <- na.omit(Control)
Kabuki <- na.omit(Kabuki)

## Filtering all the variables, that are not needed
unique(colnames(Control))

column_to_keep <- c( "Sample_ID","age..years..ch1"     ,    "disease.state.ch1"     ,  "gender.ch1" ,            
                      "sample.type.ch1"   ,      "tissue.ch1"    ,         "Basename" ,              
                      "Sentrix_ID"       ,       "Sentrix_Position"   ,"platform_id"    )
states <- list(Control, Kabuki)
states_variable_filtered <- lapply(states, function(state) state[, column_to_keep])


Controls <- states_variable_filtered[[1]]
Validation <- Control
rm(Control, Kabuki , states, column_to_keep)

m_Kabuki <-states_variable_filtered[[2]]

rm(states_variable_filtered)

m_Kabuki$Group <- m_Kabuki$disease.state.ch1
m_Kabuki$Group <- recode(m_Kabuki$Group,
                       "KDM6A variant"="Kabuki",
                       "KMT2D variant"="Kabuki")

m_Kabuki$age..years..ch1 <-as.numeric(m_Kabuki$age..years..ch1)
unique(m_Kabuki$sample.type.ch1)

Kabuki <- m_Kabuki
summary(Kabuki)

barplot(table(Kabuki$Group))

barplot(table(Kabuki$gender.ch1,Kabuki$Group))

boxplot(Kabuki$age..years..ch1 ~ Kabuki$Group)

basedir <- "/path/to/GSE97362"
file_name <- "Kabuki.csv"
full_path <- file.path(basedir, file_name)
write.csv(Kabuki,file =full_path,row.names = FALSE)

##Taking the validation cohort to be my transformed/not transformed group

Validation$age..years..ch1 <- as.numeric(Validation$age..years..ch1)
rm(Controls, file_name, n, full_path, m_Kabuki)

# selection of the group at random
set.seed(97362)
Validation$Group <- sample(c("Val1","Val2"), nrow(Validation), replace=TRUE)
Validation1 <- subset(Validation,Group == "Val1")
Validation2 <- subset(Validation,Group== "Val2")

set.seed(280704)
Validation1$Group <- sample(c("Treatment","Control"), nrow(Validation1), replace=TRUE)
Validation2$Group <- sample(c("Treatment","Control"), nrow(Validation2), replace=TRUE)

# observation of the repartition in the groups

barplot(table(Validation1$Group))
barplot(table(Validation1$gender.ch1,Validation1$Group))
boxplot(Validation1$age..years..ch1 ~ Validation1$Group)

barplot(table(Validation2$Group))
barplot(table(Validation2$gender.ch1,Validation2$Group))
boxplot(Validation2$age..years..ch1 ~ Validation2$Group)

## saving the files
file_name <- "simulation.csv"
full_path <- file.path(basedir, file_name)
write.csv(Validation2,file =full_path,row.names = FALSE)

