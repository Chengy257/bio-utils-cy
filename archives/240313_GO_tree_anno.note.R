library(plyr)
library(stringr)
library(ape)
library(GOSemSim)
library(ggtree)   ## 进化树
library(scales)
library(cowplot) 
library(ggplot2)

ggplot(df.bar,aes(x=V2,y=V3,fill=V4))+geom_bar(stat="identity")+coord_flip()+facet_wrap(~V1,nrow = 2)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")
df.bar_allgo <- read.table("../231017_GO_mainPicture/GO_BP_multiGSEA/df.gsea_gobp")
View(df.bar_allgo)
ggplot(df.bar_allgo,aes(x=V2,y=V3,fill=V4))+geom_bar(stat="identity")+coord_flip()+facet_wrap(~V1,nrow = 1)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")
pairwise_termsim(df.bar_allgo$V2)
pairwise_termsim(df.bar_allgo)
pairwise_termsim(df.bar_allgo$V2)
GOSemSim::clusterSim(df.bar_allgo$V2)
df.bar_allgo <- read.table("../231017_GO_mainPicture/GO_BP_multiGSEA/df.gsea_gobp2")
ggplot(df.bar_allgo,aes(x=V2,y=V4,fill=V5))+geom_bar(stat="identity")+coord_flip()+facet_wrap(~V3,nrow = 1)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")
ggplot(df.bar_allgo,aes(x=V2,y=V4,fill=V5))+geom_bar(stat="identity")+scale_x_discrete(label= V1)+coord_flip()+facet_wrap(~V3,nrow = 1)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")
ggplot(df.bar_allgo,aes(x=V2,y=V4,fill=V5))+geom_bar(stat="identity")+scale_x_discrete(label= df.bar_allgo$V1)+coord_flip()+facet_wrap(~V3,nrow = 1)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")
ggplot(df.bar_allgo,aes(x=V2,y=V4,fill=V5))+geom_bar(stat="identity")+coord_flip()+facet_wrap(~V3,nrow = 1)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")
library(GOSemSim)
mgoSim(unique(df.bar_allgo$V2),unique(df.bar_allgo$V2))
mgoSim(unique(df.bar_allgo$V2),unique(df.bar_allgo$V2),godata(org.Osativa.eg.db,ont = "BP"))
mgoSim(unique(df.bar_allgo$V2),unique(df.bar_allgo$V2),godata(org.Osativa.eg.db,ont = "BP"),measure = "Wang")
godata(org.Osativa.eg.db,ont = "BP")
godata(org.Osativa.eg.db,ont = "BP",keytype = "GID")
mgoSim(unique(df.bar_allgo$V2),unique(df.bar_allgo$V2),semData = mmgo,measure = "Wang")
mgoSim(unique(df.bar_allgo$V2),unique(df.bar_allgo$V2),semData = mmgo,measure = "Wang",combine = NULL)
tree <- nj(as.dist(1-ego.sim))
tree <- nj(as.dist(1-gosim))
gosim[1:3,1:3]



mmgo <- godata(org.Osativa.eg.db,ont = "BP",keytype = "GID")
gosim <- mgoSim(unique(df.bar_allgo$V2),unique(df.bar_allgo$V2),semData = mmgo,measure = "Wang",combine = NULL)

tree <- nj(as.dist(1-gosim))
p_tree <- ggtree(tree) + geom_tiplab() + #写GO term
geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + #写node编号
coord_cartesian(xlim=c(-.1,1.3)) #左右两侧留出合适的空间
barplot <- ggplot(df.bar_allgo,aes(x=V2,y=V4,fill=V5))+geom_bar(stat="identity")+coord_flip()+facet_wrap(~V3,nrow = 1)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")


library(aplot)
barplot %>% insert_left(p_tree)

p_tree <- ggtree(tree,layout = "ellipse") + geom_tiplab() + #写GO term
+     geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + #写node编号
+     coord_cartesian(xlim=c(-.1,1.3))

p_tree <- ggtree(tree,layout = "ellipse")

p_tree <- ggtree(tree,branch.length = "none")

barplot %>% insert_left(p_tree,width = .2)


ggplot(df.bar_allgo,aes(x=V2,y=V4,fill=V5))+geom_bar(stat="identity")+coord_flip()+facet_wrap(~V3,nrow = 1)+scale_fill_gradient(low = "red",high = "blue")+labs(y="NES",x=NULL,fill="p.value")


p_text <- ggplot(df.bar_allgo,aes(x=V2,y=0,label=V1))+geom_text(aes(hjust=0,family="sans",fontface ="italic"))+ylim(0,1)+theme_void()+coord_flip()

barplot %>% insert_left(p_tree,width = .2) %>% insert_right(p_text,width = .8)