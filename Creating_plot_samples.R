{library(ggplot2)
  library(plotly)
  library('gmodels')
  library(scales)
  library(reshape2)
  library(viridis)}
{library(dplyr)
  library(purrr)
  library(stringr)}

## Plot the samples gender and ages based on the groups for comparisons 


KMT2D <- read.csv( "/Users/Asus/Documents/GSEs/GSE97362/KMT2D.csv")
CHD7 <- read.csv("/Users/Asus/Documents/GSEs/GSE97362/CHD7.csv")
Simulation <- read.csv("/Users/Asus/Documents/GSEs/GSE97362/simulation.csv")

table(KMT2D$Group)
table(Simulation$Group)
table(CHD7$Group)
KMT2D_Sex <- as.data.frame(table(KMT2D$gender.ch1, KMT2D$Group))
CHD7_Sex <- as.data.frame(table(CHD7$gender.ch1,CHD7$Group))
simu_Sex <- as.data.frame(table(Simulation$gender.ch1,Simulation$Group))

simu_Sex$Var2 <- recode(simu_Sex$Var2, 
                        "control"="Control",
                        "treatment"="Treatment")

colors <- magma(2, alpha = 1, begin = 0.5, end = 0.7, direction = 1)
CHD7_Sex$Var2 <- factor(CHD7_Sex$Var2, levels = c("Control", "CHARGE"))
simu_Sex$Var2 <- factor(simu_Sex$Var2, levels = c("Control", "Treatment"))
myblanktheme <- theme(
  plot.title = element_text(family = "Arial", face = "bold", size = (15)),
  legend.title = element_text(colour = "black", face = "bold.italic", family = "Arial"),
  legend.text = element_text(face = "italic", colour = "black", family = "Arial"),
  axis.title = element_text(family = "Arial", size = (10), colour = "black"),
  axis.text =  element_text(family = "Arial", colour = "black", size = (10))
)

KMT2D_S <- ggplot(KMT2D_Sex, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = colors)
CHD7_S <- ggplot(CHD7_Sex, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = colors)
simu_S <- ggplot(simu_Sex, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = colors)

print(KMT2D_S  + myblanktheme + labs(y = "Number of observations", x="Sex", fill="Groups"))
print(simu_S  + myblanktheme + labs(y = "Number of observations", x="Sex", fill="Groups"))
print(CHD7_S  + myblanktheme + labs(y = "Number of observations", x="Sex", fill="Groups"))


CHD7$Group <- factor(CHD7$Group, levels = c("Control", "CHARGE"))

Simulation$Group <- recode(Simulation$Group, 
                        "control"="Control",
                        "treatment"="Treatment")

Simulation$Group<- factor(Simulation$Group, levels = c("Control", "Treatment"))
# plot KMT2D
gg_k <- ggplot(KMT2D, aes(x=age..years..ch1, fill= Group)) +
  geom_histogram(binwidth = 1, alpha=.5, position = "identity") +
  xlim(-1, 21) +
  ylim(0, 6)+
  scale_fill_manual(values = colors)+labs(fill = "Groups")


# plot CHD7
gg_c <- ggplot(CHD7, aes(x=age..years..ch1, fill= Group)) +
  geom_histogram(binwidth = 1, alpha=.5, position = "identity") +
  xlim(-1, 30) +
  ylim(0, 6)+
  scale_fill_manual(values = colors)+labs(fill = "Groups")


# Plot simu
gg_s <-  ggplot(Simulation, aes(x=age..years..ch1, fill= Group)) +
  geom_histogram(binwidth = 1, alpha=.5, position = "identity") +
  scale_x_continuous(limits=c(0,20)) +
  scale_y_continuous(limits=c(0,6)) +
  scale_fill_manual(values = colors) + labs(fill = "Groups")


print(gg_k  + myblanktheme + labs(y = "Number of observation", x="Age (years)"))
print(gg_c  + myblanktheme + labs(y = "Number of observation", x="Age (years)"))
print(gg_s  + myblanktheme + labs(y = "Number of observation", x="Age (years)"))

