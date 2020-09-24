# MSCWJ-1 cell movement analysis
# Danila Bobkov, 2019

# Load libraries
library(ggplot2)
library(ggsignif)

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# prepare the data
all.h$probe <- as.factor(all.h$probe)

# Set X - axis names
CellSciGuylabs <- c("9", "15", "28", "36", "38")

# Plot results
# 1 Mean speed
ggplot(all.h, aes(x = probe, y = mean_speed)) +
  
  ggtitle("Mean speed") +
  
  ylim(c(0, 120)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.2),
              cex = .5,
            shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Micrometers per hour',
       x = "Cell passage") +

  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) + 
  
#  stat_summary(fun.data=data_summary, 
#               geom="crossbar", 
#               width=0.8) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 0) +

  geom_signif(y_position = c(101),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +

  geom_signif(y_position = c(92),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04) + theme_classic(base_size=14) + 
  
  geom_signif(y_position = c(120),
              xmin = c(1),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(111),
              xmin = c(2),
              xmax = c(5),
              annotation = "ns",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(85),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(78),
              xmin = c(4),
              xmax = c(5),
              annotation = "****",
              tip_length = 0.04)


# 2 Max speed
ggplot(all.h, aes(x = probe, y = max_speed)) +
  
  ggtitle("Max speed") +
  
  ylim(c(0, 650)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.2),
              cex = .5,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Micrometers per hour',
       x = "Cell passage") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) + 
  
#  stat_summary(fun.data=data_summary, 
#               geom="crossbar", 
#               width=0.8) +

  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 0) +
  
  geom_signif(y_position = c(470),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) +

  geom_signif(y_position = c(450),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +

  geom_signif(y_position = c(630),
              xmin = c(1),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) +

  geom_signif(y_position = c(510),
              xmin = c(4),
              xmax = c(5),
              annotation = "ns",
              tip_length = 0.04) +

  geom_signif(y_position = c(570),
              xmin = c(2),
              xmax = c(4),
              annotation = "ns",
              tip_length = 0.04) +

  geom_signif(y_position = c(510),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04)
  
# 3 sinuosity

ggplot(all.h, aes(x = probe, y = sinuosity)) +
  
  ggtitle("Sinuosity") +
  
  ylim(c(0.1, 1.1)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.2),
              cex = .5,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'Sinuosity index',
       x = "Cell passage") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) + 
  
#  stat_summary(fun.data=data_summary, 
#               geom="crossbar", 
#               width=0.8) +

  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 0) +
  
  geom_signif(y_position = c(0.9),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04) +

  geom_signif(y_position = c(.83),
              xmin = c(2),
              xmax = c(3),
              annotation = "****",
              tip_length = 0.04)+
  
  geom_signif(y_position = c(0.9),
              xmin = c(3),
              xmax = c(4),
              annotation = "ns",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(1.1),
              xmin = c(1),
              xmax = c(4),
              annotation = "ns",
              tip_length = 0.04)+
  
  geom_signif(y_position = c(0.83),
              xmin = c(4),
              xmax = c(5),
              annotation = "ns",
              tip_length = 0.04)+
  
  geom_signif(y_position = c(1),
              xmin = c(1),
              xmax = c(3),
              annotation = "*",
              tip_length = 0.04)

# 4 Emax


ggplot(all.h, aes(x = probe, y = emax)) +
  
  ggtitle("Straightness") +
  
  ylim(c(0, 20)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.2),
              cex = .5,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'E-max',
       x = "Cell passage") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) + 
  
#  stat_summary(fun.data=data_summary, 
#               geom="crossbar", 
#               width=0.8) +
  
  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 0) +

  geom_signif(y_position = c(15),
              xmin = c(1),
              xmax = c(2),
              annotation = "ns",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(13),
              xmin = c(2),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(17.5),
              xmin = c(2),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(15.5),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) +
  
  geom_signif(y_position = c(19.4),
              xmin = c(4),
              xmax = c(5),
              annotation = "**",
              tip_length = 0.04)

# 5 DC

ggplot(all.h, aes(x = probe, y = DC)) +
  
  ggtitle("Directional change") +
  
  ylim(c(0.007, .048)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.2),
              cex = .5,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'DC',
       x = "Cell passage") +
  
  
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = .5, vjust = .5, 
                                   face = "plain"),
        axis.text.y = element_text(color = "black", 
                                   size = 12, angle = 0, 
                                   hjust = 1, vjust = 0, 
                                   face = "plain"),  
        axis.title.x = element_text(color = "black", 
                                    size = 12, angle = 0,
                                    hjust = .5, vjust = 0, 
                                    face = "plain"),
        axis.title.y = element_text(color = "black", 
                                    size = 12, angle = 90, 
                                    hjust = .5, vjust = .5, 
                                    face = "plain")) + 
  
#  stat_summary(fun.data=data_summary, 
#               geom="crossbar", 
#               width=0.8) +

  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 0) +
  
  geom_signif(y_position = c(.045),
              xmin = c(1),
              xmax = c(2),
              annotation = "*",
              tip_length = 0.04) + 
  
#  geom_signif(y_position = c(.049),
#              xmin = c(1),
#              xmax = c(4),
#              annotation = "****",
#              tip_length = 0.04) + 
  
  geom_signif(y_position = c(.042),
              xmin = c(2),
              xmax = c(3),
              annotation = "ns",
              tip_length = 0.04) + 
  
  geom_signif(y_position = c(.045),
              xmin = c(3),
              xmax = c(4),
              annotation = "****",
              tip_length = 0.04) + 
  
  geom_signif(y_position = c(.042),
              xmin = c(4),
              xmax = c(5),
              annotation = "***",
              tip_length = 0.04)

####### unused ggplot elements
p + geom_violin(trim = FALSE, alpha = I(0.3))


labs(title = "Mean speed", caption = "K-W p < 2.2e-16") +


# Add mean and standard deviation

p + stat_summary(fun.data=data_summary, 
             geom="crossbar", width=0.7)


p + stat_summary(fun.data=data_summary, color="blue")


