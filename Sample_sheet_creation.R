##########
# Charging all the csv created for each dataset. 

GSE97362 <- read.csv( "/Users/Asus/Documents/GSEs/GSE97362/97362_sample_sheet.csv")

############
# Load the dplyr package

{library(dplyr)
library(purrr)
library(stringr)}

# Function to separate the data into disease states

Sep <- split(GSE97362, GSE97362$sample.type.ch1)

## Grouping by type of samples

CHD7 <- rbind(Sep$`CHD7 LOF discovery cohort`, Sep$`Control for CHD7 LOF discovery cohort`)
KMT2D <- rbind(Sep$`KMT2D LOF discovery cohort`,Sep$`Control for KMT2D LOF discovery cohort`)
Control <- Sep$`Control for validation cohort`

rm(Sep, GSE97362)

# Making sure that age and Sex are present 

Control <- na.omit(Control)
CHD7 <-na.omit(CHD7)
KMT2D <- na.omit(KMT2D)

## Filtering all the variables, that are not needed
unique(colnames(Control))

column_to_keep <- c( "Sample_ID","age..years..ch1"     ,    "disease.state.ch1"     ,  "gender.ch1" ,            
                      "sample.type.ch1"   ,      "tissue.ch1"    ,         "Basename" ,              
                      "Sentrix_ID"       ,       "Sentrix_Position"   ,"platform_id"    )
states <- list(Control, CHD7, KMT2D)
states_variable_filtered <- lapply(states, function(state) state[, column_to_keep])


Controls <- states_variable_filtered[[1]]
Validation <- Control
rm(Control, CHD7, KMT2D, states, column_to_keep)



m_CHD7 <-  states_variable_filtered[[2]]
m_KMT2D <-states_variable_filtered[[3]]

rm(states_variable_filtered)

m_CHD7$Group <- m_CHD7$disease.state.ch1
m_CHD7$Group <- recode(m_CHD7$Group,
                       "CHARGE"="CHARGE",
                       "CHD7 variant"="CHARGE")
m_KMT2D$Group <- m_KMT2D$disease.state.ch1
m_KMT2D$Group <- recode(m_KMT2D$Group,
                       "KDM6A variant"="Kabuki",
                       "KMT2D variant"="Kabuki")

m_CHD7$age..years..ch1<- recode(m_CHD7$age..years..ch1, 
                              "<0.1" ="0.1")

m_CHD7$age..years..ch1<- as.numeric(m_CHD7$age..years..ch1)
m_KMT2D$age..years..ch1 <-as.numeric(m_KMT2D$age..years..ch1)
unique(m_CHD7$sample.type.ch1)
unique(m_KMT2D$sample.type.ch1)

CHD7 <- m_CHD7
KMT2D <- m_KMT2D
summary(KMT2D)
summary(CHD7)
barplot(table(KMT2D$Group))
barplot(table(CHD7$Group))
barplot(table(KMT2D$gender.ch1,KMT2D$Group))
barplot(table(CHD7$gender.ch1,CHD7$Group))
boxplot(KMT2D$age..years..ch1 ~ KMT2D$Group)
boxplot(CHD7$age..years..ch1~ CHD7$Group)

basedir <- "/Users/Asus/Documents/GSEs/GSE97362"
file_name <- "KMT2D.csv"
full_path <- file.path(basedir, file_name)
write.csv(KMT2D,file =full_path,row.names = FALSE)

file_name <- "CHD7.csv"
full_path <- file.path(basedir, file_name)
write.csv(CHD7,file =full_path,row.names = FALSE)


##Taking the validation cohort to be my transformed/not transformed group

Validation$age..years..ch1 <- as.numeric(Validation$age..years..ch1)

rm(Controls,m_CHD7, file_name, n, full_path, m_KMT2D)



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

