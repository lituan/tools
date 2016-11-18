library(reshape)
library(lattice)
sm <- read.table("similarity.txt", head = T)
md <- melt(sm, id = c("id"))
md$group[md$id <= 30] <- "Now"
md$group[md$id > 30] <- "Ancestor"
md$group[md$id == 59] <- "LCA"
md$group <- factor(md$group, levels = c("Now", "Ancestor", "LCA"))



pvalue <- c()
i <- 1
for (l in levels(md$variable)) {
  m1 <- md[md$variable == l & md$id <= 30, ]$value
  m2 <- md[md$variable == l & md$id > 40, ]$value
  h <- t.test(m1, m2,alternative = 'less')
  pvalue[i] <- h$p.value
  i <- i + 1
}
pd <- as.data.frame(list(pvalue = pvalue))
pd$indicator[pvalue <= 0.05 & pvalue > 0.01] <- "*"
pd$indicator[pvalue <= 0.01 & pvalue > 0.001] <- "**"
pd$indicator[pvalue <= 0.001 & pvalue > 1e-04] <- "***"
pd$indicator[pvalue <= 1e-04] <- "***"
pd$indicator[pvalue > 0.05] <- " "

ay <- c()
ax <- c()
i <- 1
for (l in levels(md$variable)) {
  ay[i] <- max(md[md$variable == l, ]$value) + 0.03
  ax[i] <- i
  i <- i + 1
}
pd$x <- ax
pd$y <- ay

mdd <- md[md$id <= 30 | md$id > 40,]
# stripplot(value~variable,data=md,groups=group,col=c('green','blue','red'),xlab='Repeat
# pairs',ylab='similarity',main='Distribution of Repeat
# Similarities',legend=c('now','ance','LCA'))
p <- ggplot(mdd, aes(x = variable, y = value, color = group)) + 
  geom_jitter(position = position_jitter(0.2))
p + scale_color_manual(values = c("blue", "green", "red")) + 
  annotate("text", label = pd$indicator, x = pd$x, y = pd$y) + 
  xlab("Repeat Pairs") + ylab("Similarity") + ggtitle("Distribution of Repeat similarities") + 
  theme(plot.title = element_text(hjust=0.5))

