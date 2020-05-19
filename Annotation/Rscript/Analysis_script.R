library(ggplot2)
library(ggpubr)
convert_type <- function(ms) {
  m<-ms[ms$type!="Unassigned",]
  m$Category[m$size>=50 & m$size<100] = "[50-100]"
  m$Category[m$size>=100 & m$size<250] = "]100-250]"
  m$Category[m$size>=250 & m$size<500] = "]250-500]"
  m$Category[m$size>=500 & m$size<1000] = "]500-1000]"
  m$Category[m$size>=1000] = ">1000"
  m$Category_mh[m$mh_size==0] = "0"
  m$Category_mh[m$mh_size>0 & m$mh_size<5] = "]0-5]"
  m$Category_mh[m$mh_size>5 & m$mh_size<=10] = "]5-10]"
  m$Category_mh[m$mh_size>10 & m$mh_size<=20] = "]10-20]"
  m$Category_mh[m$mh_size>20 & m$mh_size<=50] = "]20-50]"
  m$Category_mh[m$mh_size>50 & m$mh_size<=100] = "]50-100]"
  m$Category_mh[m$mh_size>100] = ">100"
  
  m$Category_loc_s[m$Repeat_loc=="SimpleRep" | m$Repeat_loc=="Simple"]= "SimpleRep"
  m$Category_loc_s[m$Repeat_loc=="Satellite"] = "SimpleRep"
  m$Category_loc_s[m$Repeat_loc=="LINE"] = "LINE"
  m$Category_loc_s[m$Repeat_loc=="SINE"] = "SINE"
  m$Category_loc_s[m$Repeat_loc=="Non repeat"] = "Other"
  m$Category_loc_s[m$Repeat_loc=="Retroposon"]= "Other TE"
  m$Category_loc_s[m$Repeat_loc=="DNA"]= "Other TE"
  m$Category_loc_s[m$Repeat_loc=="LTR"]= "Other TE"
  m$Category_loc_s[m$Repeat_loc=="RNA"] = "Other TE"
  m$Category_loc_s[m$Repeat_loc=="SegmDup"] = "SegmDup"
  #m$Category_loc_s[m$Repeat_loc=="Exon"] = "Exon"
  #m$Category_loc_s[m$Repeat_loc=="Intron"] = "Intron"
  #m$Category_loc_s[m$Repeat_loc=="Intergenic"] = "Intergenic"
  #novel, ME, TR, TD, DD, segD, Unass
  m$Category_type[m$type=="Dispersed duplication"] = "DispersDup"
  m$Category_type[m$type=="Novel sequence"] = "Novel"
  m$Category_type[m$type=="Tandem duplication"] = "TandemDup"
  m$Category_type[m$type=="Segmental duplication"] = "SegmDup"
  m$Category_type[m$type=="Tandem repeat"] = "TandemRep"
  m$Category_type[m$type=="Mobile element"] = "ME"
  m$Category_type[m$type=="Unassigned"] = "NA"
  return(m)
}

m_1 <- read.csv("~/Desktop/note_book/All_info_NA19240.csv", header = T, sep = ",",skipNul = TRUE)
otp="NA19240"

m_2 <- read.csv("~/Desktop/note_book/All_info_HG002.csv", header = T, sep = ",",skipNul = TRUE)
m_2$mh_size<-as.numeric(levels(m_2$mh_size))[m_2$mh_size] 
m_2$size<-as.numeric(levels(m_2$size))[m_2$size] 
m_2<-m_2[complete.cases(m_2), ]
otp="HG002"

m_3 <- read.csv("~/Desktop/note_book/All_info_HG002_pass.csv", header = T, sep = ",",skipNul = TRUE)
otp="HG002_pass"

m_4 <- read.csv("~/Desktop/note_book/All_info_HG00514.csv", header = T, sep = ",",skipNul = TRUE)
otp="HG00514"

m_5<- read.csv("~/Desktop/note_book/All_info_HG00733.csv", header = T, sep = ",",skipNul = TRUE)
otp="HG00733"


m_1$Individu="NA19240"
m_2$Individu="HG002"
m_3$Individu="HG002 PASS"
m_4$Individu="HG00514"
m_5$Individu="HG00733"

m_1<-convert_type(m_1)
m_2<-convert_type(m_2)
m_3<-convert_type(m_3)
m_4<-convert_type(m_4)
m_5<-convert_type(m_5)

m<-rbind(m_1,m_2,m_3,m_4,m_5)
otp="CAT_Individu"

m$Category_mh<-factor(m$Category_mh,levels=c("0","]0-10]","]10-20]","]20-50]","]50-100]",">100"))
m$Category<-factor(m$Category,levels=c("[50-100]","]100-250]","]250-500]","]500-1000]",">1000"))
m$Category_type<-factor(m$Category_type,levels=c("Novel","ME","TandemRep","TandemDup","DispersDup","NA"))
m$Category_loc_s<-factor(m$Category_loc_s,levels=c("SimpleRep","LINE","SINE","Other TE","SegmDup","Other"))
m$Individu<-factor(m$Individu,level=c("NA19240","HG00514","HG00733","HG002","HG002 PASS"))
m<-m[complete.cases(m), ]

#######################GENERAL###############################################################
g_1<-ggplot(data=m, aes(m$Category,fill=m$Individu)) + geom_histogram(stat="count", position="dodge")+ 
  labs(fill="",y="Count", x = "Insertion size (bp)") +theme(legend.position = "right",
                                                              legend.key.size = unit(2, "cm"),
                                                              legend.title = element_text(size = 25, vjust = .5, face = "bold"),
                                                              legend.text = element_text(size = 25, vjust = .5),
                                                              plot.title = element_text(hjust = 0.5, size = 25, face = "bold.italic"),
                                                              plot.subtitle = element_text(hjust = 0.5, size = 25),
                                                              plot.caption = element_text(hjust = 0, size = 25, face = "bold"),
                                                              axis.title = element_text(size = 25, face = "bold"),
                                                              axis.text.y = element_text(size = 25, face = "bold", color = "black"),
                                                              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 25, face = "bold",  color = "black"),
                                                              strip.text = element_text(size = 25, colour = "black", angle = 90, face = "bold"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) + scale_color_grey()+scale_fill_grey()

g_2<-ggplot(data=m, aes(m$Category_loc_s,fill=m$Individu)) + geom_histogram(stat="count", position="dodge")+ 
  labs(fill="",y="Count", x = "Location") +theme(legend.position = "right",
                                                 legend.key.size = unit(2, "cm"),
                                                 legend.title = element_text(size = 25, vjust = .5, face = "bold"),
                                                 legend.text = element_text(size = 25, vjust = .5),
                                                 plot.title = element_text(hjust = 0.5, size = 25, face = "bold.italic"),
                                                 plot.subtitle = element_text(hjust = 0.5, size = 25),
                                                 plot.caption = element_text(hjust = 0, size = 25, face = "bold"),
                                                 axis.title = element_text(size = 25, face = "bold"),
                                                 axis.text.y = element_text(size = 25, face = "bold", color = "black"),
                                                 axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 25, face = "bold",  color = "black"),
                                                 strip.text = element_text(size = 25, colour = "black", angle = 90, face = "bold"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) + scale_color_grey()+scale_fill_grey()

g_3<-ggplot(data=m, aes(m$Category_type,fill=m$Individu)) + geom_histogram(stat="count", position="dodge")+ 
  labs(fill="",y="Count", x = "Insertion type") +theme(legend.position = "right",
                                                       legend.key.size = unit(2, "cm"),
                                                       legend.title = element_text(size = 25, vjust = .5, face = "bold"),
                                                       legend.text = element_text(size = 25, vjust = .5),
                                                       plot.title = element_text(hjust = 0.5, size = 25, face = "bold.italic"),
                                                       plot.subtitle = element_text(hjust = 0.5, size = 25),
                                                       plot.caption = element_text(hjust = 0, size = 25, face = "bold"),
                                                       axis.title = element_text(size = 25, face = "bold"),
                                                       axis.text.y = element_text(size = 25, face = "bold", color = "black"),
                                                       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 25, face = "bold",  color = "black"),
                                                       strip.text = element_text(size = 25, colour = "black", angle = 90, face = "bold"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) + scale_color_grey()+scale_fill_grey()

g_4<-ggplot(data=m, aes(m$Category_mh,fill=m$Individu)) + geom_histogram(stat="count", position="dodge")+ 
  labs(fill="",y="Count", x = "Homology size (bp)") +theme(legend.position = "right",
                                                           legend.key.size = unit(2, "cm"),
                                                           legend.title = element_text(size = 25, vjust = .5, face = "bold"),
                                                           legend.text = element_text(size = 25, vjust = .5),
                                                           plot.title = element_text(hjust = 0.5, size = 25, face = "bold.italic"),
                                                           plot.subtitle = element_text(hjust = 0.5, size = 25),
                                                           plot.caption = element_text(hjust = 0, size = 25, face = "bold"),
                                                           axis.title = element_text(size = 25, face = "bold"),
                                                           axis.text.y = element_text(size = 25, face = "bold", color = "black"),
                                                           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 25, face = "bold",  color = "black"),
                                                           strip.text = element_text(size = 25, colour = "black", angle = 90, face = "bold"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) + scale_color_grey()+scale_fill_grey()
figure <- ggarrange(g_1, g_3, g_2,g_4, labels = c("A", "B", "C","D"),ncol = 2, nrow = 2,font.label=list(size = 25, color = "black"),common.legend = TRUE, legend = "bottom")
figure
ggsave(paste(otp,"3.png",sep=""),width = 500,
       height = 400,
       units = c("mm"),
       dpi = 600)
########################################TECHNO#####################################################

m_size <- read.csv("/home/wesley/Pub_percent_SR_size", header = T, sep = "\t",skipNul = TRUE)
m_type <- read.csv("/home/wesley/Pub_percent_SR_type_ins_noSD", header = T, sep = "\t",skipNul = TRUE)
m_loc <- read.csv("/home/wesley/Pub_percent_SR_type", header = T, sep = "\t",skipNul = TRUE)
m_mh <- read.csv("/home/wesley/Pub_percent_SR_mh", header = T, sep = "\t",skipNul = TRUE)


m_size<- subset(m_size, Individu=="NA19240" | Individu=="HG002")
m_type<- subset(m_type, Individu=="NA19240" | Individu=="HG002")
m_loc<- subset(m_loc, Individu=="NA19240" | Individu=="HG002")
m_mh<- subset(m_mh, Individu=="NA19240" | Individu=="HG002")

m_size<-m_size[complete.cases(m_size), ]
m_type<-m_type[complete.cases(m_type), ]
m_loc<-m_loc[complete.cases(m_loc), ]
m_mh<-m_mh[complete.cases(m_mh), ]

m_size$Categorys[m_size$Category=="50-100"] = "[50-100]"
m_size$Categorys[m_size$Category=="100-250"]= "]100-250]"
m_size$Categorys[m_size$Category=="250-500"]= "]250-500]"
m_size$Categorys[m_size$Category=="500-1000"] = "]500-1000]"
m_size$Categorys[m_size$Category==">1000"]= ">1000"

m_mh$Category_mhs[m_mh$Category_mh==0] = "0"
m_mh$Category_mhs[m_mh$Category_mh=="0-10"] = "]0-10]"
m_mh$Category_mhs[m_mh$Category_mh=="10-20"] = "]10-20]"
m_mh$Category_mhs[m_mh$Category_mh=="20-50"] = "]20-50]"
m_mh$Category_mhs[m_mh$Category_mh=="50-100"] = "]50-100]"
m_mh$Category_mhs[m_mh$Category_mh==">100"] = ">100"

m_size$percent[m_size$Techno=="Other technologies"]=100
m_type$percent[m_type$Techno=="Other technologies"]=100
m_loc$percent[m_loc$Techno=="Other technologies"]=100
m_mh$percent[m_mh$Techno=="Other technologies"]=100

m_size$Categorys<-factor(m_size$Categorys,levels=c("[50-100]","]100-250]","]250-500]","]500-1000]",">1000"))
m_type$Category_type<-factor(m_type$Category_type,levels=c("Novel","ME","TandemRep","TandemDup","DispersDup","SegmDup"))
m_loc$Category_loc<-factor(m_loc$Category_loc,levels=c("SimpleRep","LINE","SINE","Other TE","SegmDup","Other"))
m_mh$Category_mhs<-factor(m_mh$Category_mhs,levels=c("0","]0-10]","]10-20]","]20-50]","]50-100]",">100"))
m_size$Individu<-factor(m_size$Individu,levels=c("NA19240","HG00514","HG00733","HG002"))
m_type$Individu<-factor(m_type$Individu,levels=c("NA19240","HG00514","HG00733","HG002"))
m_loc$Individu<-factor(m_loc$Individu,levels=c("NA19240","HG00514","HG00733","HG002"))
m_mh$Individu<-factor(m_mh$Individu,levels=c("NA19240","HG00514","HG00733","HG002"))
otp="Illumina_other_compar"


SR_1<-ggplot(m_size, aes(x=m_size$Categorys, y=m_size$percent, group=m_size$Individu,fill=m_size$Individu, alpha=m_size$Techno)) + 
  geom_bar(stat="identity",position="dodge", colour="black")+ 
  labs(fill="",y="Percent", x = "Insertion size (bp)")+theme(legend.position = "right",
                                                                   legend.key.size = unit(1, "cm"),
                                                                   legend.title = element_text(size = 13, vjust = .5, face = "bold"),
                                                                   legend.text = element_text(size = 13, vjust = .5),
                                                                   plot.title = element_text(hjust = 0.5, size = 13, face = "bold.italic"),
                                                                   plot.subtitle = element_text(hjust = 0.5, size = 13),
                                                                   plot.caption = element_text(hjust = 0, size = 13, face = "bold"),
                                                                   axis.title = element_text(size = 13, face = "bold"),
                                                                   axis.text.y = element_text(size = 13, face = "bold", color = "black"),
                                                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, face = "bold",  color = "black"),
                                                                   strip.text = element_text(size = 13, colour = "black", angle = 90, face = "bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ labs(fill = "Individual", alpha="Sequencing technology") +scale_y_continuous(expand=c(0,0),limits=c(0,100))
SR_2<-ggplot(m_type, aes(x=m_type$Category_type, y=m_type$percent,group=m_type$Individu, fill=m_type$Individu, alpha=m_type$Techno)) + 
  geom_bar(stat="identity",position="dodge", colour="black")+ 
  labs(fill="",y="Percent", x = "Insertion type")+theme(legend.position = "right",
                                                             legend.key.size = unit(1, "cm"),
                                                             legend.title = element_text(size = 13, vjust = .5, face = "bold"),
                                                             legend.text = element_text(size = 13, vjust = .5),
                                                             plot.title = element_text(hjust = 0.5, size = 13, face = "bold.italic"),
                                                             plot.subtitle = element_text(hjust = 0.5, size = 13),
                                                             plot.caption = element_text(hjust = 0, size = 13, face = "bold"),
                                                             axis.title = element_text(size = 13, face = "bold"),
                                                             axis.text.y = element_text(size = 13, face = "bold", color = "black"),
                                                             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, face = "bold",  color = "black"),
                                                             strip.text = element_text(size = 13, colour = "black", angle = 90, face = "bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ labs(fill = "Individual", alpha="Sequencing technology")+scale_y_continuous(expand=c(0,0),limits=c(0,100))
SR_3<-ggplot(m_loc, aes(x=m_loc$Category_loc, y=m_loc$percent, fill=m_loc$Individu,group=m_loc$Individu, alpha=m_loc$Techno)) + 
  geom_bar(stat="identity",position="dodge", colour="black")+ 
  labs(fill="",y="Percent", x = "Insertion location")+theme(legend.position = "right",
                                                             legend.key.size = unit(1, "cm"),
                                                             legend.title = element_text(size = 13, vjust = .5, face = "bold"),
                                                             legend.text = element_text(size = 13, vjust = .5),
                                                             plot.title = element_text(hjust = 0.5, size = 13, face = "bold.italic"),
                                                             plot.subtitle = element_text(hjust = 0.5, size = 13),
                                                             plot.caption = element_text(hjust = 0, size = 13, face = "bold"),
                                                             axis.title = element_text(size = 13, face = "bold"),
                                                             axis.text.y = element_text(size = 13, face = "bold", color = "black"),
                                                             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, face = "bold",  color = "black"),
                                                             strip.text = element_text(size = 13, colour = "black", angle = 90, face = "bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ labs(fill = "Individual", alpha="Sequencing technology")+scale_y_continuous(expand=c(0,0),limits=c(0,100))
SR_4<-ggplot(m_mh, aes(x=m_mh$Category_mhs, y=m_mh$percent, fill=m_mh$Individu,group=m_mh$Individu, alpha=m_mh$Techno)) + 
  geom_bar(stat="identity",position="dodge", colour="black")+ 
  labs(fill="",y="Percent", x = "Homology size (bp)")+theme(legend.position = "right",
                                                             legend.key.size = unit(1, "cm"),
                                                             legend.title = element_text(size = 13, vjust = .5, face = "bold"),
                                                             legend.text = element_text(size = 13, vjust = .5),
                                                             plot.title = element_text(hjust = 0.5, size = 13, face = "bold.italic"),
                                                             plot.subtitle = element_text(hjust = 0.5, size = 13),
                                                             plot.caption = element_text(hjust = 0, size = 13, face = "bold"),
                                                             axis.title = element_text(size = 13, face = "bold"),
                                                             axis.text.y = element_text(size = 13, face = "bold", color = "black"),
                                                             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, face = "bold",  color = "black"),
                                                             strip.text = element_text(size = 13, colour = "black", angle = 90, face = "bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ labs(fill = "Individual", alpha="Sequencing technology")+scale_y_continuous(expand=c(0,0),limits=c(0,100))
figure <- ggarrange(SR_1, SR_2, SR_3,SR_4, labels = c("A", "B", "C","D"),ncol = 2, nrow = 2,font.label=list(size = 8, color = "black"),common.legend = TRUE, legend = "right")
figure
ggsave(paste(otp,"_comp3_all.png",sep=""),width = 400,
       height = 300,
       units = c("mm"),
       dpi = 600)


##############################################TYPE########################################################
m_1<- subset(m_1, type=!"Unassigned")
m<-m_1
m_p <- read.csv("~/Desktop/note_book/Test_Percent", header = T, sep = "\t",skipNul = TRUE)
m_p <- read.csv("~/Desktop/note_book/Formated_random.csv", header = T, sep = ",",skipNul = TRUE)


m_p<-convert_type(m_p)

m_percent<-rbind(m,m_p)
m<-m_percent
m$Category_mh<-factor(m$Category_mh,levels=c("0","]0-5]","]5-10]","]10-20]","]20-50]","]50-100]",">100"))
m$Category<-factor(m$Category,levels=c("[50-100]","]100-250]","]250-500]","]500-1000]",">1000"))
summary(m$type[m$Category_loc_s=="SegmDup"])
m$Category_type<-factor(m$Category_type,levels=c("Novel","ME","TandemRep","TandemDup","DispersDup","SegmDup","Genome percent","Random sequence"))
m$Category_loc_s<-factor(m$Category_loc_s,levels=c("SimpleRep","LINE","SINE","Other TE","SegmDup","Other"))
m$Individu<-factor(m$Individu,level=c("NA19240","HG00514","HG00733","HG002"))
f1<-ggplot(subset(m, !is.na(Category)), aes(x = Category_type, fill =Category)) + geom_bar(stat = "count",position="fill")+ 
  labs(fill="",y="Percent", x = "Insertion type") +theme(legend.position = "right",
                                                              legend.key.size = unit(2, "cm"),
                                                              legend.title = element_text(size = 25, vjust = .5, face = "bold"),
                                                              legend.text = element_text(size = 25, vjust = .5),
                                                              plot.title = element_text(hjust = 0.5, size = 25, face = "bold.italic"),
                                                              plot.subtitle = element_text(hjust = 0.5, size = 25),
                                                              plot.caption = element_text(hjust = 0, size = 25, face = "bold"),
                                                              axis.title = element_text(size = 25, face = "bold"),
                                                              axis.text.y = element_text(size = 25, face = "bold", color = "black"),
                                                              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 25, face = "bold",  color = "black"),
                                                              strip.text = element_text(size = 25, colour = "black", angle = 90, face = "bold"))+  scale_y_continuous(labels=scales::percent) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ labs(fill = "Insertion size (bp)")

f2<-ggplot(subset(m, !is.na(Category_mh)),aes(x =Category_type, fill = Category_mh)) + geom_bar(stat = "count",position="fill")+ 
  labs(fill="",y="Percent", x = "Insertion type") +theme(legend.position = "right",
                                                             legend.key.size = unit(2, "cm"),
                                                             legend.title = element_text(size = 25, vjust = .5, face = "bold"),
                                                             legend.text = element_text(size = 25, vjust = .5),
                                                             plot.title = element_text(hjust = 0.5, size = 25, face = "bold.italic"),
                                                             plot.subtitle = element_text(hjust = 0.5, size = 25),
                                                             plot.caption = element_text(hjust = 0, size = 25, face = "bold"),
                                                             axis.title = element_text(size = 25, face = "bold"),
                                                             axis.text.y = element_text(size = 25, face = "bold", color = "black"),
                                                             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 25, face = "bold",  color = "black"),
                                                             strip.text = element_text(size = 25, colour = "black", angle = 90, face = "bold"))+ scale_y_continuous(labels=scales::percent) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ labs(fill = "Homology size (bp)")
f3<-ggplot(subset(m, !is.na(Category_loc_s)), aes(x =Category_type , fill = Category_loc_s)) + geom_bar(stat = "count",position="fill")+ 
  labs(fill="",y="Percent", x = "Insertion type") +theme(legend.position = "right",
                                                   legend.key.size = unit(2, "cm"),
                                                   legend.title = element_text(size = 25, vjust = .5, face = "bold"),
                                                   legend.text = element_text(size = 25, vjust = .5),
                                                   plot.title = element_text(hjust = 0.5, size = 25, face = "bold.italic"),
                                                   plot.subtitle = element_text(hjust = 0.5, size = 25),
                                                   plot.caption = element_text(hjust = 0, size = 25, face = "bold"),
                                                   axis.title = element_text(size = 25, face = "bold"),
                                                   axis.text.y = element_text(size = 25, face = "bold", color = "black"),
                                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 25, face = "bold",  color = "black"),
                                                   strip.text = element_text(size = 25, colour = "black", angle = 90, face = "bold"))+ scale_y_continuous(labels=scales::percent)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+ labs(fill = "Location")

figure <- ggarrange(f1,f3,f2, labels = c("A", "B", "C"),ncol = 3, nrow = 1,font.label=list(size = 25, color = "black"))#,common.legend = TRUE, legend = "bottom")
figure <- ggarrange(f2,ncol = 1, nrow = 1,font.label=list(size = 25, color = "black"))#,common.legend = TRUE, legend = "bottom")

figure
ggsave(paste(otp,"Type_3.png",sep="_"),width = 300, limitsize = FALSE,
       height = 300,
       units = c("mm"),
       dpi = 600)

