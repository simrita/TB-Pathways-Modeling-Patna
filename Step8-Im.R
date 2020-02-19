##############################
# Information Matrices
##############################
iter=t
##Missing Information Matrices
CoVariance_c2=function(i,k,kd,e,f){
  return(CoVariance_2(i,k,kd,iter,e,f))
}
CoVariance_c3=function(i,k,kd,e,f,g){
  return(CoVariance_3(i,k,kd,iter,e,f,g))
}
CoVariance_c4=function(i,k,kd, e,f,g,h){
  return(CoVariance_4(i,k,kd,iter,e,f,g,h))
}

Covariance   = array(0,dim=c(n_obs,Number_Parameters_Estimated, Number_Parameters_Estimated))
Covariance_t = array(0,dim=c(n_obs, Number_Parameters_Estimated, Number_Parameters_Estimated, Stages_Interval))
Covariance_te = array(0,dim=c(n_obs, Number_Parameters_Estimated, Number_Parameters_Estimated, Stages_Interval,Stages_Interval))

I_M_t  = array(0,dim=c(n_obs,Number_Parameters_Estimated, Number_Parameters_Estimated))
I_M_te = array(0,dim=c(n_obs, Number_Parameters_Estimated, Number_Parameters_Estimated, Stages_Interval))
I_M_tem = array(0,dim=c(n_obs, Number_Parameters_Estimated, Number_Parameters_Estimated, Stages_Interval, Stages_Interval))

c=n_obs
m=Category
###### for rows corresponding to the first mix probabilities 
for(i in 1:c)
{
  
 for(e in 1:Constt0)
  {
    q=e
    k=1
    for(f in 1:Constt0)
    {
      r=f
      kd=1
      I_M_t[i,q,r]=(1/first_mix[e])*(1/first_mix[f])*CoVariance_c2(i,k,kd,e,f)-(1/first_mix[e])*(1/first_mix[m])*CoVariance_c2(i,k,kd,e,m)-(1/first_mix[m])*(1/first_mix[f])*CoVariance_c2(i,k,kd,m,f)+(1/first_mix[m])*(1/first_mix[m])*CoVariance_c2(i,k,kd,m,m)
    }
 

  


    
    q=e
    k=1
    
    for(f in 1:Category)
    {
      sum=0
      r=Constt0+f
      for(kd in 1:Stagescount(i)){
        sum=sum+(1/first_mix[e])*((Diag[i,kd]/mu_d[f])-Dur[i,kd])*CoVariance_c2(i,k,kd,e,f)-(1/first_mix[m])*((Diag[i,kd]/mu_d[f])-Dur[i,kd])*CoVariance_c2(i,k,kd,m,f)
      }
      
      I_M_t[i,q,r]=sum
    }
  

    
    q=e
    k=1
    for(f in 1:Category)
    {
      sum=0
      r=Constt1+f
      for(kd in 1:Stagescount(i)){
        sum=sum+(1/first_mix[e])*(((1-Diag[i,kd])/mu_s[f])-Dur[i,kd])*CoVariance_c2(i,k,kd,e,f)-(1/first_mix[m])*(((1-Diag[i,kd])/mu_s[f])-Dur[i,kd])*CoVariance_c2(i,k,kd,m,f)
      }
      
      I_M_t[i,q,r]=sum
    }
    
  
    q=e
    k=1
    for(f in 1:Category)
    {
      sum=0
      r=Constt2+f
      for(kd in 1:Stagescount(i)){
        sum=sum+(1/first_mix[e])*((C_Diag[i,kd]/Sens[f])-((1-C_Diag[i,kd])/(1-Sens[f])))*Diag[i,kd]*CoVariance_c2(i,k,kd,e,f)-(1/first_mix[m])*((C_Diag[i,kd]/Sens[f])-((1-C_Diag[i,kd])/(1-Sens[f])))*Diag[i,kd]*CoVariance_c2(i,k,kd,m,f)
      }
      
      I_M_t[i,q,r]=sum
    }
 


    
    q=e
    k=1
    if(Stagescount(i)>1)
    {
      for(f in 1:Category){
        
        for(g in 1:(Category-1)){
      {
        sum=0
        r=Constt3+((f-1)*(Category-1))+g
        for(kd in 1:(Stagescount(i)-1)){
          sum=sum+(1/first_mix[e])*((1-C_Diag[i,kd])/TPM[f,g])*CoVariance_c3(i,k,kd, e,f,g)-(1/first_mix[e])*((1-C_Diag[i,kd])/TPM[f,m])*CoVariance_c3(i,k,kd, e,f,m)-(1/first_mix[m])*((1-C_Diag[i,kd])/TPM[f,g])*CoVariance_c3(i,k,kd, m,f,g)+(1/first_mix[m])*((1-C_Diag[i,kd])/TPM[f,m])*CoVariance_c3(i,k,kd, m,f,m)
        }
        I_M_t[i,q,r]=sum
      }
    }
    
   
  
}
}
}
}

####### for rows corresponding to mu_d

for(i in 1:c)
{
  
  
  for(e in 1:Category)
  {
    q=Constt0+e
    for(f in 1:Constt0)
    {
      
      r=f
      kd=1
      I_M_t[i,q,r]=I_M_t[i,r,q]
    }
  

   
  
   
   
      q=Constt0+e
      for(f in 1:Category)
      {
        r=Constt0+f
        sum=0
      
        for(k in 1:Stagescount(i)) 
        {
        
        for(kd in 1:Stagescount(i)) 
         {
        sum=sum+((Diag[i,k]/(mu_d[e]))-Dur[i,k])*((Diag[i,kd]/(mu_d[f]))-Dur[i,kd])*CoVariance_c2(i,k,kd,e,f)
          }
        }
      
           I_M_t[i,q,r]=sum
          }
  }
}
  for(i in 1:c){
    for(e in 1:Category)
    {
    q=Constt0+e
    for(f in 1:Category)
    {
      r=Constt1+f
      sum=0
      
      for(k in 1:Stagescount(i)) 
      {
        
        for(kd in 1:Stagescount(i)) 
        {
          sum=sum+((Diag[i,k]/(mu_d[e]))-Dur[i,k])*(((1-Diag[i,kd])/(mu_s[f]))-Dur[i,kd])*CoVariance_c2(i,k,kd,e,f)
        }
      }
      
      I_M_t[i,q,r]=sum
    }
    }
  }
    for(i in 1:c)
    {
      for(e in 1:Category)
      {
    q=Constt0+e
    for(f in 1:Category)
    {
      r=Constt2+f
      sum=0
      
      for(k in 1:Stagescount(i)) 
      {
        
        for(kd in 1:Stagescount(i)) 
        {
          sum=sum+((Diag[i,k]/(mu_d[e]))-Dur[i,k])*((C_Diag[i,kd]/Sens[f])-(Diag[i,kd]*(1-C_Diag[i,kd])/(1-Sens[f])))*CoVariance_c2(i,k,kd,e,f)
        }
      }
      
      I_M_t[i,q,r]=sum
    }
  }
}
for(i in 1:c)
{
  for(e in 1:Category)
  {
    if(Stagescount(i)>1)
    {
    q=Constt0+e
    for(f in 1:Category)
    {
      for(g in 1:(Category-1))
      {
        
      r=Constt3+(f-1)*(Category-1)+g
      sum=0
      
      for(k in 1:(Stagescount(i)) )
      {
        
        for(kd in 1:(Stagescount(i)-1)) 
        {
          sum=sum+((Diag[i,k]/(mu_d[e]))-Dur[i,k])*((1-C_Diag[i,kd])/TPM[f,g])*CoVariance_c3(i,k,kd,e,f,g)-((Diag[i,k]/(mu_d[e]))-Dur[i,k])*((1-C_Diag[i,kd])/TPM[f,m])*CoVariance_c3(i,k,kd,e,f,m)
        }
      }
      
      I_M_t[i,q,r]=sum
    }
    }
    
  }
  }
}


# Time to Switching

for(i in 1:c)
{
  
  for(e in 1:Category)
  {
    q=Constt1+e
    for(f in 1:Constt0)
    {
      r=f
      kd=1
      I_M_t[i,q,r]=I_M_t[i,r,q]
    }
    
    
    for(f in 1:Category)
    {
      r=Constt0+f
      I_M_t[i,q,r]=I_M_t[i,r,q]
    }
    
    
   for(f in 1:Category)
   {
     r=Constt1+f
     sum=0
     for(k in 1:Stagescount(i)) 
     {
       
       for(kd in 1:Stagescount(i)) 
       {
         sum=sum+(((1-Diag[i,k])/(mu_s[e]))-Dur[i,k])*(((1-Diag[i,kd])/(mu_s[f]))-Dur[i,kd])*CoVariance_c2(i,k,kd,e,f)
       }
     }
     
     I_M_t[i,q,r]=sum
   }
   for(f in 1:Category)
   {
     r=Constt2+f
     sum=0
     
     for(k in 1:Stagescount(i)) 
     {
       
       for(kd in 1:Stagescount(i)) 
       {
         sum=sum+(((1-Diag[i,k])/(mu_s[e]))-Dur[i,k])*((C_Diag[i,kd]/Sens[f])-(Diag[i,kd]*(1-C_Diag[i,kd])/(1-Sens[f])))*CoVariance_c2(i,k,kd,e,f)
       }
     }
     
     I_M_t[i,q,r]=sum
   }
   if(Stagescount(i)>1)
   {
     q=Constt1+e
     for(f in 1:Category)
     {
       for(g in 1:(Category-1))
     {
       r=Constt3+(f-1)*(Category-1)+g
       sum=0
       
       for(k in 1:(Stagescount(i))) 
       {
         
         for(kd in 1:(Stagescount(i)-1)) 
         {
           sum=sum+(((1-Diag[i,k])/(mu_s[e]))-Dur[i,k])*((1-C_Diag[i,kd])/TPM[f,g])*CoVariance_c3(i,k,kd,e,f,g)-(((1-Diag[i,k])/(mu_s[e]))-Dur[i,k])*((1-C_Diag[i,kd])/TPM[f,m])*CoVariance_c3(i,k,kd,e,f,m)
         }
       }
       
       I_M_t[i,q,r]=sum
     }
     }
   }
  }
}
    
 
# For rows corresponding to Sensitivity

for(i in 1:c)
{
  for(e in 1:Category)
  {
    q=Constt2+e
    for(f in 1:Constt2)
    {
      r=f
      I_M_t[i,q,r]=I_M_t[i,r,q]
    }
    for(f in 1:Category)
    {
      r=Constt2+f
      sum=0
      for(k in 1:Stagescount(i)) 
      {
        
        for(kd in 1:Stagescount(i)) 
        {
          sum=sum+((C_Diag[i,k]/Sens[e])-(Diag[i,k]*(1-C_Diag[i,k])/(1-Sens[e])))*((C_Diag[i,kd]/Sens[f])-(Diag[i,kd]*(1-C_Diag[i,kd])/(1-Sens[f])))*CoVariance_c2(i,k,kd,e,f)
        }
      }
      
      I_M_t[i,q,r]=sum
    }
    if(Stagescount(i)>1)
    {
      q=Constt2+e
      for(f in 1:Category)
      {
        for(g in 1:(Category-1))
      {
        r=Constt3+(f-1)*(Category-1)+g
        sum=0
        
        for(k in 1:Stagescount(i)) 
        {
          
          for(kd in 1:(Stagescount(i)-1) )
          {
            sum=sum+((C_Diag[i,k]/Sens[e])-(Diag[i,k]*(1-C_Diag[i,k])/(1-Sens[e])))*((1-C_Diag[i,kd])/TPM[f,g])*CoVariance_c3(i,k,kd,e,f,g)-((C_Diag[i,k]/Sens[e])-(Diag[i,k]*(1-C_Diag[i,k])/(1-Sens[e])))*((1-C_Diag[i,kd])/TPM[f,m])*CoVariance_c3(i,k,kd,e,f,m)
          }
        }
        
        I_M_t[i,q,r]=sum
      }
      }
      
    }
  }
}


#for rows corresponding to TPM

for(i in 1:c)
{
  for(e in 1:Category)
  {
    for(f in 1:(Category-1))
    {
    for(g in 1:Constt3)
    {
      q=Constt3+(e-1)*(Category-1)+f
      r=g
      I_M_t[i,q,r]=I_M_t[i,r,q]
    }
    
  }
}
}


for(i in 1:c)
{
  if(Stagescount(i) > 1)
  {
  for(e in 1:Category)
  {
    for(f in 1:(Category-1))
    {
      for(g in 1:Category)
      {
        for(h in 1:(Category-1))
          {
          sum=0
          q=Constt3+(e-1)*(Category-1)+f
          r=Constt3+(g-1)*(Category-1)+h
          for(k in 1:(Stagescount(i)-1))
          {
            for(kd in 1:(Stagescount(i)-1))
            {
              sum=sum+((1-C_Diag[i,k])/TPM[e,f])*((1-C_Diag[i,kd])/TPM[g,h])*CoVariance_c4(i,k,kd, e,f,g,h)-((1-C_Diag[i,k])/TPM[e,f])*((1-C_Diag[i,kd])/TPM[g,m])*CoVariance_c4(i,k,kd,e,f,g,m)-((1-C_Diag[i,k])/TPM[e,m])*((1-C_Diag[i,kd])/TPM[g,h])*CoVariance_c4(i,k,kd, e,m,g,h)+((1-C_Diag[i,k])/TPM[e,m])*((1-C_Diag[i,kd])/TPM[g,m])*CoVariance_c4(i,k,kd, e,m,g,m)
            }
          }
          I_M_t[i,q,r]=sum
        }
      }
    }
  }
}
}



I_M=array(dim=c(Number_Parameters_Estimated, Number_Parameters_Estimated))

for(q in 1:Number_Parameters_Estimated)for(r in 1:Number_Parameters_Estimated)
{
  I_M[q, r]=sum(Observed(I_M_t[,q,r]))
}

write.csv(I_M,file="parameters/I_M.csv")























