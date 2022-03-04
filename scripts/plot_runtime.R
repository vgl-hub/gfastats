setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)

df<-read.csv("data.txt", header = TRUE, sep = "\t")

png(file="Fig 1c.png",
    width=2000, height=1000)

ggplot(df, aes(x=size, y=time, group=format)) +
  geom_point(aes(color=format), size = 3)+
  scale_color_grey() + theme_classic() +
  geom_smooth(aes(color=format)) +
  theme(
    text = element_text(size = 60),
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA),
    legend.key.size = unit(3,"cm")
  )+
  xlab("Genome size") + ylab("Time (s)")

dev.off()
