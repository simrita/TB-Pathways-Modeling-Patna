# Step-4 E-M ALgorithm



# number of iterations
iterations=100
patients=n_obs


######################################
# Step: Parameter estimation
######################################

# Defining Arrays
first_mix_A_it=array(dim=iterations)

TPM_aa_it=array(dim=iterations)
TPM_ac_it=array(dim=iterations)
TPM_ad_it=array(dim=iterations)
TPM_ba_it=array(dim=iterations)
TPM_bc_it=array(dim=iterations)
TPM_bd_it=array(dim=iterations)
TPM_ca_it=array(dim=iterations)
TPM_da_it=array(dim=iterations)



mu_d_it=array(dim=c(Category,iterations))
mu_s_it=array(dim=c(Category,iterations))
Sens_it=array(dim=c(Category,iterations))
TPM_it=array(dim=c(Category,Category,iterations))

#############
#Initial Values
#############

t=1

# first_mix_A_it[t]=length(which(Cat[,1]==1))/(length(which(Cat[,1]==1|Cat[,1]==2)))
first_mix_A_it[t]=0.5


mu_d_it[,t]=c(1/8,1/10,1/20,1/14)
mu_s_it[,t]=c(1/14,1/8,1/12,1/15)
Sens_it[,t]=c(0.5,0.3,0.2,0.6)

TPM_aa_it[t]=0.25  
TPM_ac_it[t]=0.25  
TPM_ad_it[t]=0.25  
TPM_ba_it[t]=0.25  
TPM_bc_it[t]=0.25  
TPM_bd_it[t]=0.25  
TPM_ca_it[t]=0.25  
TPM_da_it[t]=0.25  



# Paramters to be determined

first_mix_it=array(dim=c(Category,iterations))

TPM=array(dim=c(Category,Category,iterations))
S_TPM_n=array(dim=c(Category,Category,Stages-1,iterations))
W_TPM_n=array(dim=c(Category,Category,Stages-1,iterations))
S_TPM_num=array(dim=c(Category,Category,iterations))
W_TPM_num=array(dim=c(Category,Category,iterations))
TPM_num=array(dim=c(Category,Category,iterations))
S_TPM_it=array(dim=c(Category,Category,iterations))
W_TPM_it=array(dim=c(Category,Category,iterations))
TPM_num_it=array(dim=c(Category,Category,iterations))

##################### 
# E-M Algorithm
#####################

# Buliding the proportion of the mix of the patients

first_mix_calc=function(t)
{
  first_mix_c=array(dim=Category)
  
  first_mix_c[3:4]=first_mix_true[3:4]
  first_mix_c[1]=(first_mix_A_it[t])*(1-(first_mix_c[3])-(first_mix_c[4]))
  first_mix_c[2]=(1-first_mix_A_it[t])*(1-(first_mix_c[3])-(first_mix_c[4]))
  
  return(first_mix_c)
}

# Building the TPM post wrong diagnosis using the estimated paramters

TPM_calc=function(t)
{
  TPM_c=array(dim=c(Category,Category))
  
  TPM_c[3:4,3:4]=TPM_true[3:4,3:4]
  
  
  TPM_c[1,3] = (TPM_ac_it[t])
  TPM_c[1,4] = (TPM_ad_it[t])
  TPM_c[1,1] = (TPM_aa_it[t])
  
  TPM_c[2,3] = (TPM_bc_it[t])
  TPM_c[2,4] = (TPM_bd_it[t])
  TPM_c[2,1] = (TPM_ba_it[t])
  
  TPM_c[3,1] = (TPM_ca_it[t])
  TPM_c[4,1] = (TPM_da_it[t]) 
  
  TPM_c[1,2] = 1-(TPM_c[1,1])-(TPM_c[1,3])-(TPM_c[1,4])
  TPM_c[2,2] = 1-(TPM_c[2,1])-(TPM_c[2,3])-(TPM_c[2,4])
  TPM_c[3,2] = 1-(TPM_c[3,1])-(TPM_c[3,3])-(TPM_c[3,4])
  TPM_c[4,2] = 1-(TPM_c[4,1])-(TPM_c[4,3])-(TPM_c[4,4])
  
  
  return(TPM_c)
}


# Building the category indicator value for each patient at each stage

Num_Calc=function(i, e, a, t)
{
  if(e==1)
  {
    Num_c=((1-N_QC[i,a])*(1-N_KC[i,a])*  (Expz[i,a,t]))+((1-N_QC[i,a])*(N_KC[i,a])*  (N_C[i,a]))
  }
  else
  {
    if(e==2)
    {
      Num_c=((1-N_QC[i,a])*(1-N_KC[i,a])*(1-Expz[i,a,t]))+((1-N_QC[i,a])*(N_KC[i,a])*(1-N_C[i,a]))
    }
    else
    {
      if(e==3)
      {
        Num_c=(N_QC[i,a])*  (N_MC[i,a])
      }
      else
      {
        if(e==4)
        {
          Num_c=(N_QC[i,a])*(1-N_MC[i,a])
        }
        else
        {
          Num_c=0
        }
      }
    }
  }
  return(Num_c)
}

Catnew=function(i,a)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    catnew=Cat[i,a]
  }
  else
  {
    catnew=(Cat[i,a])-2
  }
  return(catnew)
}

# Possibe transitions

TPM_Count=function(i,a,t)
{
  TPMC=array(dim=c(Category,Category))
  if(Cat[i,a]==3|Cat[i,a]==4)
  {
    if((Cat[i,(a+1)]==3)|(Cat[i,(a+1)]==4))
    {
      TPMC[1,1]=Joint11(i,a,t) 
      TPMC[1,2]=Joint12(i,a,t) 
      TPMC[2,1]=Joint21(i,a,t) 
      TPMC[2,2]=Joint22(i,a,t) 
    }
    else
    {
      TPMC[1,Catnew(i,(a+1))]=  Expz[i,a,t]
      TPMC[2,Catnew(i,(a+1))]=1-Expz[i,a,t]
    }
  }
  else
  {
    if((Cat[i,(a+1)]==3)|(Cat[i,(a+1)]==4))
    {
      TPMC[Catnew(i,a),1]=  Expz[i,(a+1),t]
      TPMC[Catnew(i,a),2]=1-Expz[i,(a+1),t]
    }
    else
    {
      TPMC[Catnew(i,a),Catnew(i,(a+1))]=1
    }
  
  return(TPMC)
  }
}


#################################################################

Expz=array(dim=c(n_obs,Stages,iterations))
Num=array(dim=c(n_obs,Category,Stages,iterations))
TPM_C=array(dim=c(n_obs,Category, Category,Stages,iterations))

# setting the tolerance for log likelihood function to limit the iterations

tol = 0.0001
gap = array(dim=iterations)
logL = array(dim=iterations)
logL[1]=0

########################################
###### Running the iterations 
########################################
t=1
repeat
  {
  

  for(i in 1:n_obs)
  {
    a=1
    repeat
    {
      if(Cat[i,a]==1|Cat[i,a]==2|Cat[i,a]==5|Cat[i,a]==6)
      {
        Expz[i,a,t]=0.5
      }
      else
      {
        Expz[i,a,t]=Expzee(i,a,t)
      }
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  
  # COmputing the number of transitions
  

  for(i in 1:n_obs){
    for(e in 1:Category){
      for(f in 1:Category){
         if(Stagescount(i)>1)
           {
        a=1
        repeat
        {
          if(!is.null(dim(TPM_Count(i,a,f))))
         ##if((dim(TPM_Count(i,a,t))[1] == 4 )&& (dim(TPM_Count(i,a,t))[2] == 4))
           {
           
          TPM_C[i,e,f,a,t]=TPM_Count(i,a,t)[e,f]
           }

        
          if(a==(Stagescount(i)-1)){break}
          a=a+1
        }
      }
    }
    }
    }
  # Computing the number of patients with each category at each stage
  
  for(i in 1:n_obs) for(e in 1:Category)for(a in 1:Stages)
  {
    Num[i,e,a,t]=Num_Calc(i, e, a, t)
  }
  
  #### Estimating the parameters
  
  # Initial mix
  
  a=1
  first_mix_A_it[t+1]=(sum(Observed(Num[,1,a,t])))/((sum(Observed(Num[,1,a,t])))+(sum(Observed(Num[,2,a,t]))))
  
  # Mu_d, Mu_s and Diagnostic accuracy
  
  num=vector()
  den=vector()
  
  
  for(e in 1:Category){
    for(i in 1: patients){
      num[i]=sum(Observed(Num_Calc(i,e,,t)*Diag[i,]))
      den[i]=sum(Observed(Num_Calc(i,e,,t)*Dur[i,]))         
    }
    mu_d_it[e,t+1]=sum(num)/sum(den)
  }
  for(e in 1:Category){
    for(i in 1: patients){
      num[i]=sum(Observed(Num_Calc(i,e,,t)*(1-Diag[i,])))
      den[i]=sum(Observed(Num_Calc(i,e,,t)*Dur[i,]))         
    }
    mu_s_it[e,t+1]=sum(num)/sum(den)
  }
  for(e in 1:Category){
    for(i in 1: patients){
      num[i]=sum(Observed(Num_Calc(i,e,,t)*Diag[i,]*C_Diag[i,]))
      den[i]=sum(Observed(Num_Calc(i,e,,t)*Diag[i,]))         
    }
    Sens_it[e,t+1]=sum(num)/sum(den)
  }
  
  
  
  
  
  
  
  # TPMs 
  
  if(Stages_Interval>1)
  {
    a=1
    repeat
    {
      for(e in 1:Category)for(f in 1:Category)
      {
        W_TPM_n[e,f,a,t+1]=sum(Observed(TPM_C[,e,f,a,t]*(Diag[,a])*(1-C_Diag[,a])))
        S_TPM_n[e,f,a,t+1]=sum(Observed(TPM_C[,e,f,a,t]*(1-Diag[,a])))
        
      }  
      if(a==Stages_Interval-1){break}
      a=a+1
    }
    
    for(e in 1:Category)for(f in 1:Category)
    {
      W_TPM_num[e,f,t+1]=sum(Observed(W_TPM_n[e,f,,t+1]))
      S_TPM_num[e,f,t+1]=sum(Observed(S_TPM_n[e,f,,t+1]))
      TPM_num[e,f,t+1]=W_TPM_num[e,f,t+1]+S_TPM_num[e,f,t+1]
      
    }
    for(e in 1:Category)for(f in 1:Category)
    {
      TPM[e,f,t+1]=(TPM_num[e,f,t+1])/(sum(TPM_num[e,,t+1]))
      
    }   
  }
  
  if(Stages>1)
  {
    TPM_aa_it[t+1]=(TPM[1,1,t+1])
    TPM_ac_it[t+1]=(TPM[1,3,t+1])
    TPM_ad_it[t+1]=(TPM[1,4,t+1])
    TPM_ba_it[t+1]=(TPM[2,1,t+1])
    TPM_bc_it[t+1]=(TPM[2,3,t+1])
    TPM_bd_it[t+1]=(TPM[2,4,t+1])
    TPM_ca_it[t+1]=(TPM[3,1,t+1])
    TPM_da_it[t+1]=(TPM[4,1,t+1])
  }
    
    
    # Estimates computation
    
    first_mix_it[1:4,t+1]=first_mix_calc(t+1)[1:4]
    
    for(e in 1:Category)for(f in 1:Category)
    {
      TPM_it[e,f,t+1]=TPM_calc(t+1)[e,f]
      
    }
  
  logL1=array(dim=Category)
  logL2=array(dim=c(Category,Stages))
  logL3=array(dim=c(Category,Stages))
  logL4=array(dim=c(Category,Stages))
  logL5=array(dim=c(Category,Stages))
  logL6=array(dim=c(Category,Stages))
  logL7=array(dim=c(Category,Category))
  logL8=array(dim=c(Category,Category))
  
  for(e in 1:Category)
  {
    logL1[e]=sum(Observed((Num[,e,1,t])*(log(first_mix_it[e,t+1]))))
    for(a in 1:Stages_Interval)
    {
      logL2[e,a]=sum(Observed((Num[,e,a,t])*(Diag[,a])*  (log(mu_d_it[e,t+1]))))
      logL3[e,a]=sum(Observed((Num[,e,a,t])*(1-Diag[,a])*(log(mu_s_it[e,t+1]))))
      logL4[e,a]=sum(Observed((Num[,e,a,t])*(Dur[,a])*((mu_d_it[e,t+1])+(mu_s_it[e,t+1]))))
      logL5[e,a]=sum(Observed((Num[,e,a,t])*(Diag[,a])*(C_Diag[,a])*(log(Sens_it[e,t+1]))))
      logL6[e,a]=sum(Observed((Num[,e,a,t])*(Diag[,a])*(1-C_Diag[,a])*(log(1-Sens_it[e,t+1]))))
    }
  }
  
  for(e in 1:Category)
  {
    for(f in 1:Category)
    {
      logL7[e,f]=sum(Observed((TPM_num[e,f,t+1])*(Diag)*(1-C_Diag)*(log(TPM_it[e,f,t+1]))))
      logL8[e,f]=sum(Observed((TPM_num[e,f,t+1])*(1-Diag)         *(log(TPM_it[e,f,t+1]))))
    }
  }
  
  
  logL[t+1]=sum(Observed(logL1))+sum(Observed(logL2))+sum(Observed(logL3))-sum(Observed(logL4))+
    sum(Observed(logL5))+sum(Observed(logL6))+sum(Observed(logL7))+sum(Observed(logL8))
  
  gap[t+1]=logL[t+1]-logL[t]
  
  if(abs(gap[t+1])<=tol)
  {
    break
  }
  if(t==iterations-1)
  {
    break
  }
  
  t<-t+1
  
  
}

######## writing all the parameters estimated to the file "parameters.csv"
data=data.frame(first_mix=first_mix_it[,t],Sensitivity=Sens_it[,t],Td=1/mu_d_it[,t],Ts=1/mu_s_it[,t],TPM=TPM_it[,,t])
dir.create("parameters")
write.csv(data,"parameters/parameters.csv")


    
 