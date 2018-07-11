setwd("C:/Users/zha200/Documents/Data")
temp_drug_1=dget("widedrugdata.txt") ### Change name，temp_drug->temp_drug_1
ID_D=names(temp_drug_1)
### ID_D

temp_ADE_1=dget("wide ADE.txt")
ID_ADE=names(temp_ADE_1)
### ID_ADE

trump=intersect(ID_D,ID_ADE)
id_1=match(ID_D,trump)
id_2=match(ID_ADE,trump)

temp_drug=temp_drug_1[is.na(id_1)==F] ### Change name，temp->temp_drug
temp_ADE=temp_ADE_1[is.na(id_2)==F]

ID_D=names(temp_drug)
ID_ADE=names(temp_ADE)
sum(ID_D==ID_ADE)

h=28523
temp_drug[h]
temp_ADE[h]

ADE_names=unique(as.vector(unlist(temp_ADE)))

Drug_dic=read.csv("name.csv")
dic_ID=as.character(Drug_dic[,1])
dic_names=as.character(Drug_dic[,3])

##########################↑↑↑↑↑Read Data↑↑↑↑↑############################
######################↓↓↓↓↓All Drugs + All ADE↓↓↓↓↓#######################

drug_all=rep(0,1764)
ADE_all=rep(0,15351)
drug_ADE_all=matrix(0,15351,1764)

for(i in 1:4077447)
{
  try_drug=temp_drug[[i]]
  try_drug_name=dic_names[match(try_drug,dic_ID)]
  
  try_ADE=temp_ADE[[i]]
  
  x_ADE=match(try_ADE,ADE_names)
  y_drug=match(try_drug_name,dic_names)
  ADE_all[x_ADE]=ADE_all[x_ADE]+1
  drug_all[y_drug]=drug_all[y_drug]+1
  drug_ADE_all[x_ADE,y_drug]=drug_ADE_all[x_ADE,y_drug]+1
  
  cat("iter=",i,"\n")
}

colnames(drug_ADE_all)=dic_names
write.csv(drug_ADE_all,"drug_ADE_all.csv",row.names=ADE_names)

names(drug_all)=dic_names
names(ADE_all)=ADE_names

write.csv(drug_all,"drug_all.csv")
write.csv(ADE_all,"ADE_all.csv")

sum=Reduce("+",drug_ADE_all)
write.csv(sum,"drug_ADE_sum.csv",row.names=F)

######################↓↓↓↓↓47 drug-drug pairs↓↓↓↓↓#######################

######################↓↓↓↓↓N11D & N111↓↓↓↓↓#######################

drug_1=read.csv("drug_1.csv")
drug_2=read.csv("drug_2.csv")
drug_pair=data.frame(drug_1,drug_2)
drug_pair
drug_pair_vector=paste(drug_pair[,1],drug_pair[,2])
drug_pair_vector
N11D=matrix(0,15351,47)
N111=matrix(0,15351,47)

for(j in 1:4077447)
{
  try_drug=temp_drug[[j]]  ### drug information
  try_drug_name=dic_names[match(try_drug,dic_ID)]
  
  a=rep(try_drug_name, length(try_drug_name))
  b=rep(try_drug_name, each=length(try_drug_name))
  c=paste(a,b)
  
  temp_drug_pair=match(c,drug_pair_vector)
  temp_drug_pair=temp_drug_pair[!is.na(temp_drug_pair)]
  
  try_ADE=temp_ADE[[j]] ### ADE information
  
  temp_drug_ADE=match(try_ADE,ADE_names)
  
  N11D[,temp_drug_pair]=N11D[,temp_drug_pair]+1
  N111[temp_drug_ADE,temp_drug_pair]=N111[temp_drug_ADE,temp_drug_pair]+1
  
  cat("iter =",j,"\n")
}

colnames(N11D)=drug_pair_vector
write.csv(N11D,"N11D.csv",row.names=ADE_names)

colnames(N111)=drug_pair_vector
write.csv(N111,"N111.csv",row.names=ADE_names)

######################↓↓↓↓↓N10D & N101↓↓↓↓↓#######################

N10D=matrix(0,15351,47)
N101=matrix(0,15351,47)

drug_1_vector=as.character(drug_1$x)
drug1_index_in_all=match(drug_1_vector, dic_names)

######convert drug_all from vector to matrix######
drug_all_t=t(drug_all)
drug_all_t_rep=drug_all_t[rep(1, times = 15351), ]
##################################################

N10D=drug_all_t_rep[,drug1_index_in_all]
N10D=N10D-N11D

N101=drug_ADE_all[,drug1_index_in_all]
N101=N101-N111

colnames(N10D)=drug_pair_vector
write.csv(N10D,"N10D.csv",row.names=ADE_names)
colnames(N101)=drug_pair_vector
write.csv(N101,"N101.csv",row.names=ADE_names)

######################↓↓↓↓↓N01D & N011↓↓↓↓↓#######################

N01D=matrix(0,15351,47)
N011=matrix(0,15351,47)

drug_2_vector=as.character(drug_2$y)
drug2_index_in_all=match(drug_2_vector, dic_names)

N01D=drug_all_t_rep[,drug2_index_in_all]
N01D=N01D-N11D

N011=drug_ADE_all[,drug2_index_in_all]
N011=N011-N111

colnames(N01D)=drug_pair_vector
write.csv(N01D,"N01D.csv",row.names=ADE_names)
colnames(N011)=drug_pair_vector
write.csv(N011,"N011.csv",row.names=ADE_names)

######################↓↓↓↓↓N00D & N001↓↓↓↓↓#######################

N00D=matrix(0,15351,47)
N001=matrix(0,15351,47)
all_matrix=matrix(length(ID_D),15351,47)

######convert ADE_all from vector to matrix######
ADE_all_rep=replicate(47, ADE_all)
##################################################

N00D=all_matrix-N01D-N10D-N11D
N001=ADE_all_rep-N011-N101-N111

colnames(N00D)=drug_pair_vector
write.csv(N00D,"N00D.csv",row.names=ADE_names)
colnames(N001)=drug_pair_vector
write.csv(N001,"N001.csv",row.names=ADE_names)

########################↓↓↓↓↓2 Drugs + 1 ADE↓↓↓↓↓#########################
##############################Verification###############################

firstDrug=("FLUCONAZOLE")
secondDrug=("METHADONE")
both=c("FLUCONAZOLE","METHADONE")

drug_combo_ADE=matrix(0,4,2)

ADE=("RHABDOMYOLYSIS")

for(i in 1:4077447)
{
  try_drug=temp_drug[[i]]  ### drug information
  try_drug_name=dic_names[match(try_drug,dic_ID)]
  
  temp_firstDrug=1-is.na(match(firstDrug,try_drug_name))
  temp_secondDrug=1-is.na(match(secondDrug,try_drug_name))
  temp_neither=0+(sum(is.na(match(both,try_drug_name)))==2)
  temp_both=0+(sum(is.na(match(both,try_drug_name)))==0)
  
  drug_combo=matrix(0,4,2)
  drug_combo[1,]=temp_neither
  drug_combo[2,]=temp_firstDrug-temp_both
  drug_combo[3,]=temp_secondDrug-temp_both
  drug_combo[4,]=temp_both
  
  try_ADE=temp_ADE[[i]] ### ADE information
  temp_drug_ADE=1-is.na(match(ADE,try_ADE)) #match-->match
  
  drug_combo[,1]=drug_combo[,1]*temp_drug_ADE
  
  drug_combo_ADE=drug_combo_ADE+drug_combo
  
  cat("iter=",i,"\n")
}

drug_combo_ADE

########################Omega Calculation#########################

f00=N001/N00D
f10=N101/N10D
f01=N011/N01D
f11=N111/N11D

EXP_1=pmax((f00/(1-f00)),(f01/(1-f01)))+pmax((f00/(1-f00)),(f10/(1-f10)))-(f00/(1-f00))+1
EXP_2=1-1/EXP_1
EXP=EXP_2*(N11D)

omega_0_unfiltered=log((N111/EXP),base=2)
write.csv(omega_0_unfiltered,"omega_0_unfiltered.csv",row.names=ADE_names)

omega_unfiltered=log(((N111+0.5)/(EXP+0.5)),base=2)
write.csv(omega_unfiltered,"omega_unfiltered.csv",row.names=ADE_names)

omega_sd_unfiltered=1/(N111*(log(2))^2)
omega_025_unfiltered=omega_0_unfiltered-1.96*omega_sd_unfiltered
write.csv(omega_025_unfiltered,"omega_025_unfiltered.csv",row.names=ADE_names)

#########################Cluster Analysis#########################
#################Remove -Inf/Inf/NA cols and rows#################
#############################omega_0#############################

rownames(omega_0_unfiltered)=ADE_names
omega_0_unfiltered[omega_0_unfiltered==-Inf] <- 0
omega_0_unfiltered[omega_0_unfiltered==Inf] <- 0

omega_col=rep(0,ncol(omega_0_unfiltered))
omega_row=rep(0,nrow(omega_0_unfiltered))

for(m in 1:ncol(omega_0_unfiltered))
{
  omega_col[m]=sum(omega_0_unfiltered[,m])
}

for(n in 1:nrow(omega_0_unfiltered))
{
  omega_row[n]=sum(omega_0_unfiltered[n,])
}

omega_0=omega_0_unfiltered[omega_row!=0,omega_col!=0]
dim(omega_0)

write.csv(omega_0,"omega_0.csv")

##############################omega##############################

rownames(omega_unfiltered)=ADE_names
omega_unfiltered[omega_unfiltered==-Inf] <- 0
omega_unfiltered[omega_unfiltered==Inf] <- 0

omega_col=rep(0,ncol(omega_unfiltered))
omega_row=rep(0,nrow(omega_unfiltered))

for(m in 1:ncol(omega_unfiltered))
{
  omega_col[m]=sum(omega_unfiltered[,m])
}

for(n in 1:nrow(omega_unfiltered))
{
  omega_row[n]=sum(omega_unfiltered[n,])
}

omega=omega_unfiltered[omega_row!=0,omega_col!=0]
dim(omega)

write.csv(omega,"omega.csv")

############################omega_025############################

rownames(omega_025_unfiltered)=ADE_names
omega_025_unfiltered[omega_025_unfiltered==-Inf] <- 0
omega_025_unfiltered[omega_025_unfiltered==Inf] <- 0

omega_col=rep(0,ncol(omega_025_unfiltered))
omega_row=rep(0,nrow(omega_025_unfiltered))

for(m in 1:ncol(omega_025_unfiltered))
{
  omega_col[m]=sum(omega_025_unfiltered[,m])
}

for(n in 1:nrow(omega_025_unfiltered))
{
  omega_row[n]=sum(omega_025_unfiltered[n,])
}

omega_025=omega_025_unfiltered[omega_row!=0,omega_col!=0]
dim(omega_025)

write.csv(omega_025,"omega_025.csv")

###############################Plot###############################

setwd("C:/Users/zha200/Documents/Data/Graph")
library(gplots)
library(RColorBrewer)
library(randomcoloR)

#########################test_unfiltered#########################

pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(1,15,length.out=2))
##length(pairs.breaks)
colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")

jpeg("test.jpg", width=16, height=30, units="in", res=1000)
hc <- hclust(dist(omega_0))
groups <- cutree(hc, k=100) 
rc <- brewer.pal(12, "Paired")
rc <- colorRampPalette(rc)(100)
test_col_name=colnames(omega_0[,colSums(omega_0>1)>0])
test_row_name=rownames(omega_0[rowSums(omega_0>1)>0,])
res <- heatmap.2(omega_0,breaks=pairs.breaks,col=colfunc,
                 RowSideColors=rc[groups], srtCol=45,
                 main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                 cexCol=1,cexRow=0.1,
                 margins=c(15,6),keysize=0.5,
                 labRow=test_row_name,labCol=test_col_name)
dev.off()

Order_M=omega_0[,(res$colInd)]
write.csv(Order_M,"Order_Result_test.csv")

#################Another method to assign colors#################

# pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(1, 15,length.out=2))
# ##length(pairs.breaks)
# colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
# 
# jpeg("test.jpg", width=18, height=24, units="in", res=1000)
# hc <- hclust(dist(omega_0))
# groups <- cutree(hc, k=15) 
# rc <- distinctColorPalette(15)
# res <- heatmap.2(omega_0,breaks=pairs.breaks,col=colfunc,
#                  RowSideColors=rc[groups],
#                  main="",trace="none",cexCol=1.5,cexRow=0.1,margins=c(30,30))
# dev.off()
# 
# Order_M=omega_0[,(res$colInd)]
# write.csv(Order_M,"Order_Result_test.csv")

############################|omega_0|############################
############################|omega_0|############################
############################|omega_0|############################

###############################d>1###############################

f=1
d=1
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)

  temp_omega=omega_0[(rowSums(omega_0>f))>thres_omega[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>1.5###############################

f=1
d=1.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_0[(rowSums(omega_0>f))>thres_omega[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>2###############################

f=1
d=2
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_0[(rowSums(omega_0>f))>thres_omega[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>2.5###############################

f=1
d=2.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega)) 
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_0[(rowSums(omega_0>f))>thres_omega[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>3###############################

f=1
d=3
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups],srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_0[(rowSums(omega_0>f))>thres_omega[i],]
  temp_name=paste("fig_omega_0_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_0_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups],srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

############################Combo d>1############################

f=1
d=1
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>1.5############################

f=1
d=1.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>2############################

f=1
d=2
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>2.5############################

f=1
d=2.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>3############################

f=1
d=3
setwd("/Users/Lulu/Documents/Intern/Graph/omega_0")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_0[ADE_all[rownames(omega_0)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

#############################|omega|#############################
#############################|omega|#############################
#############################|omega|#############################

###############################d>1###############################

f=1
d=1
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega[(rowSums(omega>f))>thres_omega[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>1.5###############################

f=1
d=1.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega[(rowSums(omega>f))>thres_omega[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>2###############################

f=1
d=2
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega[(rowSums(omega>f))>thres_omega[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>2.5###############################

f=1
d=2.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega[(rowSums(omega>f))>thres_omega[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>3###############################

f=1
d=3
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega[(rowSums(omega>f))>thres_omega[i],]
  temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=44) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(44)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

############################Combo d>1############################

f=1
d=1
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

###########################Combo d>1.5###########################

f=1
d=1.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>2############################

f=1
d=2
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

###########################Combo d>2.5###########################

f=1
d=2.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>3############################

f=1
d=3
setwd("/Users/Lulu/Documents/Intern/Graph/omega")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega[ADE_all[rownames(omega)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

###########################|omega_025|###########################
###########################|omega_025|###########################
###########################|omega_025|###########################

###############################d>1###############################

f=1
d=1
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_025[(rowSums(omega_025>f))>thres_omega[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>1.5###############################

f=1
d=1.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_025[(rowSums(omega_025>f))>thres_omega[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>2###############################

f=1
d=2
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_025[(rowSums(omega_025>f))>thres_omega[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>2.5###############################

f=1
d=2.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_025[(rowSums(omega_025>f))>thres_omega[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

###############################d>3###############################

f=1
d=3
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19)
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)
  
  temp_omega=omega_025[(rowSums(omega_025>f))>thres_omega[i],]
  temp_name=paste("fig_omega_025_f",f,"_d",d,"_",thres_omega[i],".jpeg",sep="")
  order_result_name=paste("Order_Result_omega_025_f",f,"_d",d,"_",thres_omega[i],".csv",sep="")
  pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
  colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
  jpeg(temp_name, width=16, height=30, units="in", res=1000)
  hc <- hclust(dist(temp_omega))
  groups <- cutree(hc, k=19) 
  rc <- brewer.pal(12, "Paired")
  rc <- colorRampPalette(rc)(19)
  test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
  test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
  res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                   RowSideColors=rc[groups], srtCol=45,
                   main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                   cexCol=1,cexRow=0.1,
                   margins=c(15,6),keysize=0.5,
                   labRow=test_row_name,labCol=test_col_name)
  dev.off()
  Order_M=temp_omega[,(res$colInd)]
  write.csv(Order_M,order_result_name)  
  
  cat("iter=",i,"\n")
}

############################Combo d>1############################

f=1
d=1
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>1.5############################

f=1
d=1.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>2############################

f=1
d=2
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>2.5############################

f=1
d=2.5
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}

############################Combo d>3############################

f=1
d=3
setwd("/Users/Lulu/Documents/Intern/Graph/omega_025")

thres_ADE_freq=c(1000,3000,5000,7000,9000,15000,20000,25000,30000,40000,50000)
thres_omega=c(0,1,2,3,4,5,6,7,8,9,10)

for(i in 1:11)
{
  temp_omega=omega_025[ADE_all[rownames(omega_025)]>thres_ADE_freq[i],]
  for(j in 1:11)
  {
    temp_omega=temp_omega[(rowSums(temp_omega>f))>thres_omega[j],]
    temp_name=paste("fig_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".jpeg",sep="")
    order_result_name=paste("Order_Result_omega_f",f,"_d",d,"_",thres_omega[j],"_ADE_freq_",thres_ADE_freq[i],".csv",sep="")
    
    pairs.breaks <- c(seq(-15, 0,length.out = 2),seq(d,15,length.out=2))
    colfunc<- colorpanel(n=3,low="white",mid="grey",high="black")
    jpeg(temp_name, width=16, height=30, units="in", res=1000)
    hc <- hclust(dist(temp_omega))
    groups <- cutree(hc, k=2)
    rc <- brewer.pal(12, "Paired")
    rc <- colorRampPalette(rc)(2)
    test_col_name=colnames(temp_omega[,colSums(temp_omega>d)>0])
    test_row_name=rownames(temp_omega[rowSums(temp_omega>d)>0,])
    res <- heatmap.2(temp_omega,breaks=pairs.breaks,col=colfunc,
                     RowSideColors=rc[groups], srtCol=45,
                     main="",xlab="Drug Pairs",ylab="ADE",trace="none",
                     cexCol=1,cexRow=0.1,
                     margins=c(15,6),keysize=0.5,
                     labRow=test_row_name,labCol=test_col_name)
    dev.off()
    Order_M=temp_omega[,(res$colInd)]
    write.csv(Order_M,order_result_name)
    
    cat("iter=",i,"-",j,"\n")
  }
}
