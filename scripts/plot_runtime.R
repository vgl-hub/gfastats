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
    legend.key.size = unit(3,"cm"),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
  ) +
  xlab("Genome size (Gbp)") + ylab("Time (s)") +
  guides(color=guide_legend(override.aes=list(fill=NA)))

dev.off()
