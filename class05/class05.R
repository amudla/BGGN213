# ...
# "title: Class 5: Graphic"

# class 5 R grapics and plots
# output: 
old.par <- par()
#secetion 1A
data <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot(x=data$Age,y=data$Weight,type = 'o', pch=15,cex=1.5,lwd=2,ylim=c(2,10),xlab="Age (months)",ylab="Weight (kg)",main="Body Weight with Age")

#section 2B

data2 <- read.delim("bimm143_05_rstats/feature_counts.txt",header = TRUE)
par(mar=c(3.1, 11.1, 4.1, 2))
barplot(data2$Count,horiz=TRUE,xlab ="Count",names.arg = data2$Feature,main="Number of Counts",las=1)
#par(mar=c(3.1, 11.1, 4.1, 2),yaxs = "r")

# section 2C
par(mar=c(5, 5, 2, 2))
data3 <- c(rnorm(10000),rnorm(10000)+4)
hist(data3,breaks = 100)

# section 3A
data4 <- read.delim("bimm143_05_rstats/male_female_counts.txt",header = TRUE)
par(mar=c(6, 5, 2, 2))
barplot(data4$Count,col = rainbow(10),names.arg = data4$Sample,las=2,ylab = "counts")


# section 3B
data5  <- read.delim("bimm143_05_rstats/up_down_expression.txt",header = TRUE)
str(data5)
table(data5$State)
par(mar=c(5, 5, 2, 2))
plot(x=data5$Condition1,y=data5$Condition2,col= data5$State,xlab = "condition1",ylab = "condition2",main = "Expression Level")
levels(data5$State)
palette(c("blue","grey","red"))

# section 3C
data6 <- read.delim("bimm143_05_rstats/expression_methylation.txt",header = TRUE)
index <- data6$expression>0
color <- c("blue","green","red")
densecol <- densCols(x=data6$gene.meth[index],y=data6$expression[index],colramp = colorRampPalette(color))
plot(x=data6$gene.meth[index],y=data6$expression[index],col = densecol,pch = 20,xlab = "Gene Methylation", ylab = "Expression",main = "Effect of Gene Methylation on Gene Expression")

# section 4A
plot(x=data6$promoter.meth,y=data6$gene.meth,xlab = "Promoter Methylation",ylab = "Gene Methylation")
source("bimm143_05_rstats/color_to_value_map.r")

colorRamp <- colorRampPalette(c("blue","red"))(100)
value <- data6$expression
PerPointcolor <- map.colors(value,c(min(value),max(value)),colorRamp)
plot(x=data6$promoter.meth,y=data6$gene.meth,xlab = "Promoter Methylation",ylab = "Gene Methylation",col = PerPointcolor)

# # Supplement 2: Color 1,2,3
# set.seed(19)
# X <- rnorm(30)
# Y <- rnorm(30)
# plot(X, Y, col = rep(1:3, each = 10), pch = 19)
# legend("bottomright", legend = paste("Group", 1:3), col = 1:3, pch = 19, bty = "n")
# 
# # playing with Color
# library(RColorBrewer)
# display.brewer.all()