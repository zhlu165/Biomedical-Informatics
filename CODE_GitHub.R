setwd("C:/Users/zha200/Documents/Data")
temp_drug_1=dget("widedrugdata.txt") ### Change name，temp_drug->temp_drug_1
ID_D=names(temp_drug_1)
### ID_D

temp_ADE_1=dget("wide ADE.txt")
ID_ADE=names(temp_ADE_1)
### ID_ADE

trump=intersect(ID_D,ID_ADE)
id_1=pmatch(ID_D,trump)
id_2=pmatch(ID_ADE,trump)

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
  try_drug_name=dic_names[pmatch(try_drug,dic_ID)]
  
  try_ADE=temp_ADE[[i]]
  
  x_ADE=pmatch(try_ADE,ADE_names)
  y_drug=pmatch(try_drug_name,dic_names)
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
  try_drug_name=dic_names[pmatch(try_drug,dic_ID)]
  
  a=rep(try_drug_name, length(try_drug_name))
  b=rep(try_drug_name, each=length(try_drug_name))
  c=paste(a,b)
  
  temp_drug_pair=pmatch(c,drug_pair_vector)
  temp_drug_pair=temp_drug_pair[!is.na(temp_drug_pair)]
  
  try_ADE=temp_ADE[[j]] ### ADE information
  
  temp_drug_ADE=pmatch(try_ADE,ADE_names)
  
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
drug1_index_in_all=pmatch(drug_1_vector, dic_names, duplicate=1)

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
drug2_index_in_all=pmatch(drug_2_vector, dic_names, duplicate=1)

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

########################Omega Calculation#########################

f00=N001/N00D
f10=N101/N10D
f01=N011/N01D
f11=N111/N11D

EXP_1=pmax((f00/(1-f00)),(f01/(1-f01)))+pmax((f00/(1-f00)),(f10/(1-f10)))-(f00/(1-f00))+1
EXP_2=1-1/EXP_1
EXP=EXP_2*(N11D)
Omega=log((N111/EXP),base=2)
write.csv(Omega,"Omega.csv")
Omega_sd=1/(N111*(log(2))^2)
write.csv(Omega_sd,"Omega_sd.csv")
Omega_025=Omega-1.96*Omega_sd
write.csv(Omega_025,"Omega_025.csv")

########################↓↓↓↓↓2 Drugs + 1 ADE↓↓↓↓↓#########################

simv=("FLUCONAZOLE")
niso=("METHADONE")
both=c("FLUCONAZOLE","METHADONE")

drug_combo_ADE=matrix(0,4,2)

ADE=("RHABDOMYOLYSIS")

for(i in 1:4077447)
{
  try_drug=temp_drug[[i]]  ### drug information
  try_drug_name=dic_names[pmatch(try_drug,dic_ID)]
  
  temp_simv=1-is.na(pmatch(simv,try_drug_name))  ### simv yes or no
  temp_niso=1-is.na(pmatch(niso,try_drug_name))
  temp_neither=0+(sum(is.na(pmatch(both,try_drug_name)))==2)
  temp_both=0+(sum(is.na(pmatch(both,try_drug_name)))==0)
  
  drug_combo=matrix(0,4,2)
  drug_combo[1,]=temp_neither
  drug_combo[2,]=temp_simv-temp_both
  drug_combo[3,]=temp_niso-temp_both
  drug_combo[4,]=temp_both
  
  try_ADE=temp_ADE[[i]] ### ADE information
  temp_drug_ADE=1-is.na(pmatch(ADE,try_ADE))
  
  drug_combo[,1]=drug_combo[,1]*temp_drug_ADE
  
  drug_combo_ADE=drug_combo_ADE+drug_combo
  
  cat("iter=",i,"\n")
}

drug_combo_ADE
