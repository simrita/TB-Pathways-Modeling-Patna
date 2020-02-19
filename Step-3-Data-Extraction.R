# This fragment of code picks the variables from the Mumbai data and creates the pathways in terms of observables and unobservables.

# Step-I Initialization

# Basic Parameters

ptm <- proc.time()

n_obs=Sample_Size_2.0_TB
Categories=6
Category=4
Stages=7
Stages_Interval=5

#Defining arrays

Dur=array(dim=c(n_obs, Stages))
Cat=array(dim=c(n_obs, Stages))
Ca=array(dim=c(n_obs, Stages))
Diag=array(dim=c(n_obs, Stages))
C_Diag=array(dim=c(n_obs, Stages))
N_C=array(dim=c(n_obs,Stages))
N_KC=array(dim=c(n_obs,Stages))
N_QC=array(dim=c(n_obs,Stages))
N_MC=array(dim=c(n_obs,Stages))
#Simulated pathway

#User-defined functions

Categor=function(i,k)
{
  if(Cat[i,k]==1|Cat[i,k]==3)
  {
    Categ=1
  }
  else
  {
    if(Cat[i,k]==2|Cat[i,k]==4)
    {
      Categ=2
    }
    else
    {
      Categ = (Cat[i,k])-2
    }
  }
  return(Categ)
}

# Indicator variable for providers with known qualifications as Public or Chemists

N_QCat=function(i,k)
{
  if(Cat[i,k]==5|Cat[i,k]==6)
  {
    N_QCate=1
  } 
  else
  {
    N_QCate=0
  } 
  return(N_QCate)
}

# Indicator variable for providers with public facility 

N_MCat=function(i,k)
{
  if(Cat[i,k]==5)
  {
    N_MCate=1
  } 
  else
  {
    N_MCate=0
  } 
  return(N_MCate)
}

# Indicator variables for providers with known qualifications as IFP and FP

N_KCat=function(i,k)
{
  if(Cat[i,k]==1|Cat[i,k]==2)
  {
    N_KCate=1
  } 
  else
  {
    N_KCate=0
  } 
  return(N_KCate)
}

# Indicator variables for providers with known qualifications as IFP

N_Cat=function(i,k)
{
  if(Cat[i,k]==1)
  {
    N_Cate=1
  } 
  else
  {
    N_Cate=0
    
  } 
  return(N_Cate)
}

# Indicator for Diagnosis

Diagn=function(i,k)
{
  if(Indicator_Provider_Diagnosis[i,k]=="Yes")
  {
    dia=1
  }
  else
  {
    dia=0
  }    
  return(dia)
}

# Indicator for Correct Diagnosis

Correct_Diag=function(i,k)
{
  if(Diagnosis_Given_By_Provider[i,k]=="TB"|Diagnosis_Given_By_Provider[i,k]=="MDR TB")
  {
      diagn=1
  }
  else
  {
    diagn=0
  }
  return(diagn)
}

# Number of Stages till first correct diagnosis for each patient 
Stagescount=function(i)
{
  k=1
  repeat
  {
    if(C_Diag[i,k]==1|k==Stages)
    {
      Stagec=k
    }
    else
    {
      Stagec=0
    }
    if(C_Diag[i,k]==1|k==Stages){break}
    k=k+1
  }
  return(Stagec)
}  

# Creation of the pathway with the indicator variables until first correct diagnosis

for(i in 1:n_obs)
{
  k=1
  repeat  
  {
    # Category
      Cat[i,k]=Categorynew[i,k]
      Ca[i,k]=Categor(i,k)
    # Duration with the provider
      Dur[i,k]=Time_Category[i,k]
    # Diagnosis by the provider
      Diag[i,k]=Diagn(i,k)
    # Correct Diagnosis by the provider
      C_Diag[i,k]=Correct_Diag(i,k)
    # Indicator variable for qualification and facility
      N_QC[i,k]=N_QCat(i,k)
      N_C[i,k]=N_Cat(i,k)
      N_KC[i,k]=N_KCat(i,k)
      N_MC[i,k]=N_MCat(i,k)

    if(k==Stages|C_Diag[i,k]==1){break}
    k=k+1
  }
}

Pathways= cbind(Category=Cat[,1],Dur=Dur[,1],Diag=Diag[,1],C_Diag=C_Diag[,1],
                Category=Cat[,2],Dur=Dur[,2],Diag=Diag[,2],C_Diag=C_Diag[,2],
                Category=Cat[,3],Dur=Dur[,3],Diag=Diag[,3],C_Diag=C_Diag[,3],
                Category=Cat[,4],Dur=Dur[,4],Diag=Diag[,4],C_Diag=C_Diag[,4],
                Category=Cat[,5],Dur=Dur[,5],Diag=Diag[,5],C_Diag=C_Diag[,5])

#######################################
# Step-3 Paramater's true values
#######################################

# Events at a Stage

first_mix_true=array(dim=Category)

Sens_true=array(dim=Category)

for(e in 1:Category)
{
  k=1
  first_mix_true[e]=length(which(Ca[,k]==e))/n_obs
  
  Sens_true[e]= length(which((Ca[,]==e)&(C_Diag[,]==1)))/
                    length(which((Ca[,]==e)&  (Diag[,]==1)))
}

# Transition across stages

TPM_true= array(dim=c(Category,Category))
TPM_n_true= array(dim=c(Category,Category))
TP_true= array(dim=c(Category,Category,Stages-1))
SP_true= array(dim=c(Category,Category,Stages-1))
WP_true= array(dim=c(Category,Category,Stages-1))
# Count of all transitions per stage

k=1
repeat
{
  for(e in 1:Category)for(f in 1:Category)
  {
    SP_true[e,f,k]=length(which(((Ca[,k]==e)&(Diag[,k]==0))&(Ca[,k+1]==f)))
    WP_true[e,f,k]=length(which(((Ca[,k]==e)&(Diag[,k]==1 & C_Diag[,k]==0))&(Ca[,k+1]==f)))
    TP_true[e,f,k]=SP_true[e,f,k]+WP_true[e,f,k]
  }
  if(k==Stages-1){break}
  k=k+1
}

# Count of transitions across stage

for(e in 1:Category)for(f in 1:Category)
{
  TPM_n_true[e,f]=sum(TP_true[e,f,])
}

# Computation of TPM

for(e in 1:Category)for(f in 1:Category)
{

  TPM_true[e,f]=TPM_n_true[e,f]/sum(TPM_n_true[e,])
  
}

# Parameter values

first_mix_true
Sens_true
TPM_true
proc.time() - ptm