##### ANALYSIS OF SBP BEFORE AND DURING AORTIC OCLUSSION USING REBOA
### IN PATIENTS WITH ACRETISMO

### Data SBP
rm(list = ls())
library(readxl)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggsci)
library(gridExtra)

database <- read_excel("data/REBOA_Acreta_Database.xlsx", 
                                    sheet = "Sheet3")


database$delta1 <- (database$pos1 - database$p1) / database$p1



#### Without hypotensive
a1 <- stack(database[5,2:29])
a2 <- stack(database[6,2:29])
a3 <- stack(database[7,2:29])
a4 <- stack(database[9,2:29])
a5 <- stack(database[10,2:29])

lancet <- pal_lancet("lanonc")(10)
miscolores <- colorRampPalette(lancet)


graph <- ggplot(data = a1)+
  geom_line(aes(x = ind, y = values, group = 1), color = lancet[1])+
  geom_point(aes(x = ind, y = values), color = lancet[1])


graph <- graph +
  geom_line(data = a2, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[2])+
  geom_point(data = a2, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[2])

graph <- graph +
  geom_line(data = a3, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[3])+
  geom_point(data = a3, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[3])

graph <- graph +
  geom_line(data = a4, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[4])+
  geom_point(data = a4, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[4])

graph <- graph +
  geom_line(data = a5, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[5])+
  geom_point(data = a5, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[5])

hypo1 <- graph+
  theme_cowplot()+
  geom_vline(xintercept = 2, linetype = 2)+
  geom_vline(xintercept = 3, linetype = 2)+
  geom_vline(xintercept = 13, linetype = 2)+
  geom_vline(xintercept = 14, linetype = 2)+
  geom_vline(xintercept = 18, linetype = 2)+
  geom_vline(xintercept = 19, linetype = 2)+
  geom_hline(yintercept = 90, linetype = 3, alpha = 0.7)+
  geom_hline(yintercept = 60, linetype = 3, alpha = 0.7)+
  ylab("Diastolic Blood Pressure - mm Hg")+
  xlab("")+
  scale_y_continuous(limits = c(40,90), breaks = seq(40,90,by = 10))+
  scale_x_discrete(labels = c("BL","T0", "T1", "T2",
                                "T3", "T4", "T5", "T6",
                                "T7", "T8", "T9", "T10",
                                "T15", "L5", "L4", "L3",
                              "L2", "L1","P1","P2","P3",
                              "P4","P5","P6","P7","P8",
                              "P9","P15") ) 
  
hypo1
#### With hypotensive
a1 <- stack(database[1,2:29])
a2 <- stack(database[2,2:29])
a3 <- stack(database[3,2:29])
a4 <- stack(database[4,2:29])
a5 <- stack(database[8,2:29])

lancet <- pal_lancet("lanonc")(9)
miscolores <- colorRampPalette(lancet)


graph <- ggplot(data = a1)+
  geom_line(aes(x = ind, y = values, group = 1), color = lancet[6])+
  geom_point(aes(x = ind, y = values), color = lancet[6])


graph <- graph +
  geom_line(data = a2, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[7])+
  geom_point(data = a2, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[7])

graph <- graph +
  geom_line(data = a3, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[8])+
  geom_point(data = a3, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[8])

graph <- graph +
  geom_line(data = a4, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[9])+
  geom_point(data = a4, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[9])

graph <- graph +
  geom_line(data = a5, aes(x = ind, y = values, group = 1), inherit.aes = FALSE, color = lancet[1])+
  geom_point(data = a5, aes(x = ind, y = values), inherit.aes = FALSE, color = lancet[1])

hypo2 <-graph+
  theme_cowplot()+
  geom_vline(xintercept = 2, linetype = 2)+
  geom_vline(xintercept = 3, linetype = 2)+
  geom_vline(xintercept = 13, linetype = 2)+
  geom_vline(xintercept = 14, linetype = 2)+
  geom_vline(xintercept = 18, linetype = 2)+
  geom_vline(xintercept = 19, linetype = 2)+
  geom_hline(yintercept = 60, linetype = 3, alpha = 0.7)+
  geom_hline(yintercept = 90, linetype = 3, alpha = 0.7)+
  ylab("Diastolic Blood Pressure - mm Hg")+
  xlab("")+
  scale_y_continuous(limits = c(40,90), breaks = seq(40,90,by = 10))+
  scale_x_discrete(labels = c("BL","T0", "T1", "T2",
                              "T3", "T4", "T5", "T6",
                              "T7", "T8", "T9", "T10",
                              "T15", "L5", "L4", "L3",
                              "L2", "L1","P1","P2","P3",
                              "P4","P5","P6","P7","P8",
                              "P9","P15") ) 

png("figures/PAD.png", res = 300)
grid.arrange(hypo1, hypo2, nrow= 1)
dev.off()
