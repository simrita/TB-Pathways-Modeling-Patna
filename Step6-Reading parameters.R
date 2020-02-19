# THis fragment of code reads the estimates usng EM algorithm already stored in a csv file

# Reading the parameter files

ptm <- proc.time()
#Category=4
#n_obs=76
#Stages_Interval=5
#first_mix_temp=read.csv("parameters/First_mix.csv")
#T_d_temp=read.csv("parameters/T_d.csv")
#T_s_temp=read.csv("parameters/T_s.csv")
#sens_temp=read.csv("parameters/Sensitivity.csv")
TPM_temp=read.csv("parameters/TPM.csv")


first_mix=array(dim=Category)

mu_d=array(dim=c(Category,Stages_Interval))
mu_s=array(dim=c(Category,Stages_Interval))
Sens=array(dim=c(Category,Stages_Interval))

TPM=array(dim=c(Category,Category))


first_mix = first_mix_temp[,-1]
mu_d      = 1/T_d_temp[,-1]
mu_s      = 1/T_s_temp[,-1]
Sens      = sens_temp[,-1]
TPM     = TPM_temp[,-1]



# Number of Patients with a provider category at a stage

Num_1_temp=read.csv("parameters/Num_1.csv")
Num_2_temp=read.csv("parameters/Num_2.csv")
Num_3_temp=read.csv("parameters/Num_3.csv")
Num_4_temp=read.csv("parameters/Num_4.csv")
Num_5_temp=read.csv("parameters/Num_5.csv")



Numb_1_temp=array(dim=c(n_obs, Category))
Numb_2_temp=array(dim=c(n_obs, Category))
Numb_3_temp=array(dim=c(n_obs, Category))
Numb_4_temp=array(dim=c(n_obs, Category))
Numb_5_temp=array(dim=c(n_obs, Category))

Numb_1_temp=Num_1_temp[,-1]
Numb_2_temp=Num_2_temp[,-1]
Numb_3_temp=Num_3_temp[,-1]

Numb_4_temp=Num_4_temp[,-1]
Numb_5_temp=Num_5_temp[,-1]

Numb=array(dim=c(n_obs, Category, Stages_Interval))



for(i in 1:n_obs)for(e in 1:Category)
{
  Numb[i,e,1]=Numb_1_temp[i,e]
  Numb[i,e,2]=Numb_2_temp[i,e]
  Numb[i,e,3]=Numb_3_temp[i,e]
  Numb[i,e,4]=Numb_4_temp[i,e]
  Numb[i,e,5]=Numb_5_temp[i,e]
}

Expz_temp=read.csv("parameters/Expz.csv")

Expz_t=Expz_temp[,-1]

Exp_z=array(dim=c(n_obs, Stages_Interval))
for(k in 1:Stages_Interval)
{
  Exp_z[,k]=Expz_t[,k]
}

# Number of transitions 

TPM_num_temp=read.csv("parameters/TPM_num.csv")

TPM_numb=TPM_num_temp[,-1]




proc.time() - ptm