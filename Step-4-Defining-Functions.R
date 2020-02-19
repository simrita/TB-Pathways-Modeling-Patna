#####################
# User-defined functions
#########################

ptm <- proc.time()

#### Probability of Td<=Ts
iterations=100
mu_d_it=array(dim=c(Category,iterations))
mu_s_it=array(dim=c(Category,iterations))

pr_diag=function(e,a,t)
{
  pr_d= (mu_d_it[e,t])/((mu_d_it[e,t])+(mu_s_it[e,t]))
  return(pr_d)
}

pr=function(i,e,a,t)
{
  pr= dexp(Dur[i,a], ((mu_d_it[e,t])+(mu_s_it[e,t])))
  return(pr)
}

##### To create a set of binaries with digits x

binary <-  function(x, b=2)
{
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N <- length(x)
  xMax <- max(x)	
  ndigits <- (floor(logb(xMax, base=2))+1)
  Base.b <- array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits)
  {#i <- 1
    Base.b[, ndigits-i+1] <- (x %% b)
    x <- (x %/% b)
  }
  if(N ==1) Base.b[1, ] else Base.b
}

##Binary possibilities corresponding to m stages

ZW=function(m)
{
  dime=2^m
  ZWa=array(dim=c(dime,m))
  for(f in 1:m)for(e in 1:dime)
  {
    ZWa[e,f]=binary(x=0:(dime-1))[e,f]
  }
  return(ZWa)
}

# Favourable permutations for ith patient

Z=function(i)
{
  Zdee=ZW(Stagescount(i))
  a=1
  repeat
  {
    if(Cat[i,a]==1)
    {
      Zde=subset(Zdee, Zdee[,a]==1)
    }
    else
    {
      if(Cat[i,a]==2)
      {
        Zde=subset(Zdee,Zdee[,a]==0)
      }
      else
      {
        if(Cat[i,a]==5)
        {
          Zde=subset(Zdee, Zdee[,a]==1)
          # Zde=replace(Zde, Cat==5, 5)
          Zde[,a]=5
        }
        else
        {
          if(Cat[i,a]==6)
          {
            Zde=subset(Zdee, Zdee[,a]==1)
            # Zde=replace(Zde,Cat==6, 6)
            Zde[,a]=6
          }
          else
          {
            Zde=Zdee
            
          }
        }
      }
    }
    Zdee=Zde
    if(a==Stagescount(i)){break}
    a=a+1
  }
  return(Zde)
}

Zd=function(i,n)
{
  ZX=subset(Z(i), (Z(i)[,n]==1))
  return(ZX)
}

# Associated category with the favourable permutations

Cat_num_sim=function(i,a,m,n)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    {
      cats=((Zd(i,n)[m,a])*1)+((1-Zd(i,n)[m,a])*2)
    }
  }
  return(cats)
}
Cat_den_sim=function(i,a,m)
{
  if(Cat[i,a]==1|Cat[i,a]==2)
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    {
      cats=((Z(i)[m,a])*1)+((1-Z(i)[m,a])*2)
    }
  }
  return(cats)
}

###############################

##################################

# favourable permutations for transitions

Hd11=function(i,n)
{
  if(Stagescount(i)>1)
  {
    Hdtemp=subset(Z(i), ((Z(i)[,n]==1)&(Z(i)[,(n+1)]==1)))
  }
  return(Hdtemp)              
}

Hd12=function(i,n)
{
  if(Stagescount(i)>1)
  {
    Hdtemp=subset(Z(i), ((Z(i)[,n]==1)&(Z(i)[,(n+1)]==0)))
  }
  return(Hdtemp)              
}

Hd21=function(i,n)
{
  if(Stagescount(i)>1)
  {
    Hdtemp=subset(Z(i), ((Z(i)[,n]==0)&(Z(i)[,(n+1)]==1)))
  }
  return(Hdtemp)              
}

Hd22=function(i,n)
{
  if(Stagescount(i)>1)
  {
    Hdtemp=subset(Z(i), ((Z(i)[,n]==0)&(Z(i)[,(n+1)]==0)))
  }
  return(Hdtemp)              
}

# Associated Categories

Cat_num_sim_h_11=function(i,a,m,n)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    { 
      cats=((Hd11(i,n)[m,a])*1)+((1-Hd11(i,n)[m,a])*2)
    }
  }
  return(cats)
}

Cat_num_sim_h_22=function(i,a,m,n)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    { 
      cats=((Hd22(i,n)[m,a])*1)+((1-Hd22(i,n)[m,a])*2)
    }
  }
  return(cats)
}

Cat_num_sim_h_12=function(i,a,m,n)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    { 
      cats=((Hd12(i,n)[m,a])*1)+((1-Hd12(i,n)[m,a])*2)
    }
  }
  return(cats)
}

Cat_num_sim_h_21=function(i,a,m,n)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    { 
      cats=((Hd21(i,n)[m,a])*1)+((1-Hd21(i,n)[m,a])*2)
    }
  }
  return(cats)
}

# Favourabl permutations for transitions for information matrix_1

Hkd11=function(i,n,nd)
{
  Hdtemp=subset(Z(i), (Z(i)[,n]==1)&(Z(i)[,nd]==1))
  return(Hdtemp)              
}

Cat_num_sim_hd_11=function(i,a,m,n,nd)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    {
      cats=((Hkd11(i,n,nd)[m,a])*1)+((1-Hkd11(i,n,nd)[m,a])*2)
    }
  }
  return(cats)
}

Cat_den_sim_h=function(i,a,m)
{
  if(Cat[i,a]==1|Cat[i,a]==2)
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    {
      cats=((Z(i)[m,a])*1)+((1-Z(i)[m,a])*2)
    }
  }
  return(cats)
}

# Favourabl permutations for transitions for information matrix_2

Hkdf11=function(i,n,nd)
{
  Hdtemp=subset(Z(i), ((Z(i)[,n]==1)&(Z(i)[,nd]==1)&(Z(i)[,(nd+1)]==1)))
  return(Hdtemp)              
}

Cat_num_sim_hdf_11=function(i,a,m,n,nd)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    {
      cats=((Hkdf11(i,n,nd)[m,a])*1)+((1-Hkdf11(i,n,nd)[m,a])*2)
    }
  }
  return(cats)
}

# Favourabl permutations for transitions for information matrix_3

Hkds11=function(i,n,nd)
{
  Hdtemp=subset(Z(i), ((Z(i)[,nd]==1|Z(i)[,nd]==0)&(Z(i)[,n]==1|Z(i)[,n]==0)&(Z(i)[,(n+1)]==1|Z(i)[,(n+1)]==0)))
  return(Hdtemp)              
}

Cat_num_sim_hds_11=function(i,a,m,n,nd)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    {
      cats=((Hkds11(i,n,nd)[m,a])*1)+((1-Hkds11(i,n,nd)[m,a])*2)
    }
  }
  return(cats)
}

# Favourabl permutations for transitions for information matrix_4

Hkdw11=function(i,n,nd)
{
  Hdtemp=subset(Z(i),((Z(i)[,n]==1)&(Z(i)[,(n+1)]==1)&
                         (Z(i)[,nd]==1)&(Z(i)[,(nd+1)]==1)))
  return(Hdtemp)              
}

Cat_num_sim_hdw_11=function(i,a,m,n,nd)
{
  if((Cat[i,a]==1)|(Cat[i,a]==2))
  {
    cats=Cat[i,a]
  }
  else
  {
    if((Cat[i,a]==5)|(Cat[i,a]==6))
    {
      cats=(Cat[i,a])-2
    }
    else
    {
      cats=((Hkdw11(i,n,nd)[m,a])*1)+((1-Hkdw11(i,n,nd)[m,a])*2)
    }
  }
  return(cats)
}

#######################
# Final building of Expz
#######################

# m is the number of rows
# n is the number of column
# t is the number of iteration

first_mix_num=function(i,m,n,t)
{
  a=1
  first_mix_n=first_mix_calc(t)[Cat_num_sim(i,a,m,n)]
  return(first_mix_n)
}

first_mix_den=function(i,m,n,t)
{
  a=1
  first_mix_d=first_mix_calc(t)[Cat_den_sim(i,a,m,n)]
  return(first_mix_d)
}

Prob_num=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim(i,a,m,n),a,t)*   pr_diag(Cat_num_sim(i,a,m,n),a,t)*  
                 (Sens_it[Cat_num_sim(i,a,m,n),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim(i,a,m,n),a,t)*   pr_diag(Cat_num_sim(i,a,m,n),a,t)*
                 (1-Sens_it[Cat_num_sim(i,a,m,n),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])*pr(i,Cat_num_sim(i,a,m,n),a,t)*(1-pr_diag(Cat_num_sim(i,a,m,n),a,t)))
  }
  return(probab1)
}

Prob_den=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab2=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_den_sim(i,a,m,n),a,t)*   pr_diag(Cat_den_sim(i,a,m,n),a,t)*  
                 (Sens_it[Cat_den_sim(i,a,m,n),t]))
    }
    else
    {
      probab2=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_den_sim(i,a,m,n),a,t)*   pr_diag(Cat_den_sim(i,a,m,n),a,t)*
                 (1-Sens_it[Cat_den_sim(i,a,m,n),t]))
    }
  }
  else
  {
    probab2=((1-Diag[i,a])                *pr(i,Cat_den_sim(i,a,m,n),a,t)*(1-pr_diag(Cat_den_sim(i,a,m,n),a,t)))
    
  }
  return(probab2)
}

Tprob_num=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim(i,a,m,n)),(Cat_num_sim(i,(a+1),m,n))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim(i,a,m,n)),(Cat_num_sim(i,(a+1),m,n))])
  }
  return(tprobab1)
}

Tprob_den=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab2= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_den_sim(i,a,m,n)),(Cat_den_sim(i,(a+1),m,n))])
    }
  }
  else
  {
    tprobab2= ((1-Diag[i,a])*TPM_calc(t)[(Cat_den_sim(i,a,m,n)),(Cat_den_sim(i,(a+1),m,n))])
  }
  return(tprobab2)
}

primer_num=function(i,m,n,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num(i,m,n,t)*Prob_num(i,a,m,n,t)
  }
  else
  {
    a=1
    num=first_mix_num(i,m,n,t)*Prob_num(i,a,m,n,t)
    a=2
    repeat
    {
      num=num * (Tprob_num(i,a-1,m,n,t)*Prob_num(i,a,m,n,t))
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer_den=function(i,m,n,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    den=first_mix_den(i,m,n,t)*Prob_den(i,a,m,n,t)
  }
  else
  {
    a=1
    den=first_mix_den(i,m,n,t)*Prob_den(i,a,m,n,t)
    a=2
    repeat
    {
      den=den*(Tprob_den(i,a-1,m,n,t)*Prob_den(i,a,m,n,t))
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=den
  return(check)
}

primer2_num=function(i,n,t)
{
  m=1
  num1=0
  repeat
  {
    num1=num1+primer_num(i,m,n,t)
    if(m==(dim(Zd(i,n))[[1]])){break}
    m=m+1
  }
  check2=num1
  return(check2)
}

primer2_den=function(i,n,t)
{
  m=1
  den1=0
  repeat
  {
    den1=den1+primer_den(i,m,n,t)
    if(m==(dim(Wd(i,n))[[1]])){break}
    m=m+1
  }
  check2=den1
  return(check2)
}

Expzee=function(i,n,t)
{
  expzee=(primer2_num(i,n,t))/((primer2_den_h(i,t)))
  return(expzee)
}

#########################################################
# Similar set of computations for transition matrices
#########################################################

first_mix_den_h=function(i,m,t)
{
  a=1
  first_mix_d=first_mix_calc(t)[Cat_den_sim_h(i,a,m)]
  return(first_mix_d)
}

Prob_den_h=function(i,a,m,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab2=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_den_sim_h(i,a,m),a,t)*   pr_diag(Cat_den_sim_h(i,a,m),a,t)* 
                 (Sens_it[Cat_den_sim_h(i,a,m),t]))
    }
    else
    {
      probab2=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_den_sim_h(i,a,m),a,t)*   pr_diag(Cat_den_sim_h(i,a,m),a,t)*
                 (1-Sens_it[Cat_den_sim_h(i,a,m),t]))
      
    }
  }
  else
  {
    probab2=((1-Diag[i,a])*pr(i,Cat_den_sim_h(i,a,m),a,t)*(1-pr_diag(Cat_den_sim_h(i,a,m),a,t)))
    
  }
  return(probab2)
}

Tprob_den_h=function(i,a,m,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab2= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_den_sim_h(i,a,m)),(Cat_den_sim_h(i,(a+1),m))])
    }
  }
  else
  {
    tprobab2= ((1-Diag[i,a])*TPM_calc(t)[(Cat_den_sim_h(i,a,m)),(Cat_den_sim_h(i,(a+1),m))])
    
  }
  return(tprobab2)
}

primer_den_h=function(i,m,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    den=first_mix_den_h(i,m,t)*Prob_den_h(i,a,m,t)
  }
  else
  {
    a=1
    den=first_mix_den_h(i,m,t)*Prob_den_h(i,a,m,t)
    a=2
    repeat
    {
      den=den*Tprob_den_h(i,(a-1),m,t)*Prob_den_h(i,a,m,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=den
  return(check)
}

primer2_den_h=function(i,t)
{
  m=1
  den1=0
  repeat
  {
    den1=den1+primer_den_h(i,m,t)
    if(m==(dim(Z(i))[[1]])){break}
    m=m+1
  }
  check2=den1
  return(check2)
}


#########################################################
# Similar set of computations for transition matrices
#########################################################

first_mix_num_h_11=function(i,m,n,t)
{
  a=1
  first_mix_n=first_mix_calc(t)[Cat_num_sim_h_11(i,a,m,n)]
  return(first_mix_n)
}

first_mix_num_h_22=function(i,m,n,t)
{
  a=1
  first_mix_n=first_mix_calc(t)[Cat_num_sim_h_22(i,a,m,n)]
  return(first_mix_n)
}

first_mix_num_h_12=function(i,m,n,t)
{
  a=1
  first_mix_n=first_mix_calc(t)[Cat_num_sim_h_12(i,a,m,n)]
  return(first_mix_n)
}

first_mix_num_h_21=function(i,m,n,t)
{
  a=1
  first_mix_n=first_mix_calc(t)[Cat_num_sim_h_21(i,a,m,n)]
  return(first_mix_n)
}

Prob_num_h_11=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_h_11(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_11(i,a,m,n),a,t)* 
                 (Sens_it[Cat_num_sim_h_11(i,a,m,n),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_h_11(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_11(i,a,m,n),a,t)*
                 (1-Sens_it[Cat_num_sim_h_11(i,a,m,n),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_h_11(i,a,m,n),a,t)*(1-pr_diag(Cat_num_sim_h_11(i,a,m,n),a,t)))
  }
  return(probab1)
}

Prob_num_h_12=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_h_12(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_12(i,a,m,n),a,t)* 
                 (Sens_it[Cat_num_sim_h_12(i,a,m,n),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_h_12(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_12(i,a,m,n),a,t)*
                 (1-Sens_it[Cat_num_sim_h_12(i,a,m,n),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_h_12(i,a,m,n),a,t)*(1-pr_diag(Cat_num_sim_h_12(i,a,m,n),a,t)))
  }
  return(probab1)
}

Prob_num_h_21=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_h_21(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_21(i,a,m,n),a,t)* 
                 (Sens_it[Cat_num_sim_h_21(i,a,m,n),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_h_21(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_21(i,a,m,n),a,t)*
                 (1-Sens_it[Cat_num_sim_h_21(i,a,m,n),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_h_21(i,a,m,n),a,t)*(1-pr_diag(Cat_num_sim_h_21(i,a,m,n),a,t)))
  }
  return(probab1)
}

Prob_num_h_22=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_h_22(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_22(i,a,m,n),a,t)* 
                 (Sens_it[Cat_num_sim_h_22(i,a,m,n),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_h_22(i,a,m,n),a,t)*   pr_diag(Cat_num_sim_h_22(i,a,m,n),a,t)*
                 (1-Sens_it[Cat_num_sim_h_22(i,a,m,n),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_h_22(i,a,m,n),a,t)*(1-pr_diag(Cat_num_sim_h_22(i,a,m,n),a,t)))
  }
  return(probab1)
}

Tprob_num_h_11=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_h_11(i,a,m,n)),(Cat_num_sim_h_11(i,(a+1),m,n))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_h_11(i,a,m,n)),(Cat_num_sim_h_11(i,(a+1),m,n))])
    
  }
  return(tprobab1)
}

Tprob_num_h_22=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_h_22(i,a,m,n)),(Cat_num_sim_h_22(i,(a+1),m,n))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_h_22(i,a,m,n)),(Cat_num_sim_h_22(i,(a+1),m,n))])
    
  }
  return(tprobab1)
}

Tprob_num_h_12=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_h_12(i,a,m,n)),(Cat_num_sim_h_12(i,(a+1),m,n))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_h_12(i,a,m,n)),(Cat_num_sim_h_12(i,(a+1),m,n))])
    
  }
  return(tprobab1)
}

Tprob_num_h_21=function(i,a,m,n,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_h_21(i,a,m,n)),(Cat_num_sim_h_21(i,(a+1),m,n))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_h_21(i,a,m,n)),(Cat_num_sim_h_21(i,(a+1),m,n))])
    
  }
  return(tprobab1)
}

primer_num_h_11=function(i,m,n,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_h_11(i,m,n,t)*Prob_num_h_11(i,a,m,n,t)
  }
  else
  {
    a=1
    num=first_mix_num_h_11(i,m,n,t)*Prob_num_h_11(i,a,m,n,t)
    a=2
    repeat
    {
      num=num*Tprob_num_h_11(i,(a-1),m,n,t)*Prob_num_h_11(i,a,m,n,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer_num_h_12=function(i,m,n,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_h_12(i,m,n,t)*Prob_num_h_12(i,a,m,n,t)
  }
  else
  {
    a=1
    num=first_mix_num_h_12(i,m,n,t)*Prob_num_h_12(i,a,m,n,t)
    a=2
    repeat
    {
      num=num*Tprob_num_h_12(i,(a-1),m,n,t)*Prob_num_h_12(i,a,m,n,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer_num_h_21=function(i,m,n,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_h_21(i,m,n,t)*Prob_num_h_21(i,a,m,n,t)
  }
  else
  {
    a=1
    num=first_mix_num_h_21(i,m,n,t)*Prob_num_h_21(i,a,m,n,t)
    a=2
    repeat
    {
      num=num*Tprob_num_h_21(i,(a-1),m,n,t)*Prob_num_h_21(i,a,m,n,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer_num_h_22=function(i,m,n,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_h_22(i,m,n,t)*Prob_num_h_22(i,a,m,n,t)
  }
  else
  {
    a=1
    num=first_mix_num_h_22(i,m,n,t)*Prob_num_h_22(i,a,m,n,t)
    a=2
    repeat
    {
      num=num*Tprob_num_h_22(i,(a-1),m,n,t)*Prob_num_h_22(i,a,m,n,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer2_num_h_11=function(i,n,t)
{
  m=1
  num1=0
  repeat
  {
    num1=num1+primer_num_h_11(i,m,n,t)
    if(m==(dim(Hd11(i,n))[[1]])){break}
    m=m+1
  }
  check2=num1
  return(check2)
}

primer2_num_h_12=function(i,n,t)
{
  m=1
  num1=0
  repeat
  {
    num1=num1+primer_num_h_12(i,m,n,t)
    if(m==(dim(Hd12(i,n))[[1]])){break}
    m=m+1
  }
  check2=num1
  return(check2)
}

primer2_num_h_21=function(i,n,t)
{
  m=1
  num1=0
  repeat
  {
    num1=num1+primer_num_h_21(i,m,n,t)
    if(m==(dim(Hd21(i,n))[[1]])){break}
    m=m+1
  }
  check2=num1
  return(check2)
}

primer2_num_h_22=function(i,n,t)
{
  m=1
  num1=0
  repeat
  {
    num1=num1+primer_num_h_22(i,m,n,t)
    if(m==(dim(Hd22(i,n))[[1]])){break}
    m=m+1
  }
  check2=num1
  return(check2)
}

Joint11=function(i,n,t)
{
  j=(primer2_num_h_11(i,n,t))/((primer2_den_h(i,t)))
  return(j)
}

Joint12=function(i,n,t)
{
  j=(primer2_num_h_12(i,n,t))/((primer2_den_h(i,t)))
  return(j)
}
Joint21=function(i,n,t)
{
  j=(primer2_num_h_21(i,n,t))/((primer2_den_h(i,t)))
  return(j)
}
Joint22=function(i,n,t)
{
  j=(primer2_num_h_22(i,n,t))/((primer2_den_h(i,t)))
  return(j)
}

#########################################################
# Similar set of computations for Information matrices
#########################################################

first_mix_num_hd_11=function(i,m,n,nd,t)
{
  a=1
  first_mix_n=first_mix_calc(t)[Cat_num_sim_hd_11(i,a,m,n,nd)]
  return(first_mix_n)
}

Prob_num_hd_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_hd_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hd_11(i,a,m,n,nd),a,t)* 
                 (Sens_it[Cat_num_sim_hd_11(i,a,m,n,nd),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_hd_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hd_11(i,a,m,n,nd),a,t)*
                 (1-Sens_it[Cat_num_sim_hd_11(i,a,m,n,nd),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_hd_11(i,a,m,n,nd),a,t)*(1-pr_diag(Cat_num_sim_hd_11(i,a,m,n,nd),a,t)))
  }
  return(probab1)
}

Tprob_num_hd_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_hd_11(i,a,m,n,nd)),(Cat_num_sim_hd_11(i,(a+1),m,n,nd))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_hd_11(i,a,m,n,nd)),(Cat_num_sim_hd_11(i,(a+1),m,n,nd))])
    
  }
  return(tprobab1)
}

primer_num_hd_11=function(i,m,n,nd,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_hd_11(i,m,n,nd,t)*Prob_num_hd_11(i,a,m,n,nd,t)
  }
  else
  {
    a=1
    num=first_mix_num_hd_11(i,m,n,nd,t)*Prob_num_hd_11(i,a,m,n,nd,t)
    a=2
    repeat
    {
      num=num*Tprob_num_hd_11(i,a-1,m,n,nd,t)*Prob_num_hd_11(i,a,m,n,nd,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer2_num_hd_11=function(i,n,nd,t)
{
  
  numh1=0
  m=1
  repeat
  {
    numh1=numh1+primer_num_hd_11(i,m,n,nd,t)
    if(m==(dim(Hkd11(i,n,nd))[[1]])){break}
    m=m+1
  }
  check2=numh1
  return(check2)
}

Joint_hd_11=function(i,n,nd,t)
{
  j=(primer2_num_hd_11(i,n,nd,t))/((primer2_den_h(i,t)))
  return(j)
}

#########################################################
# Similar set of computations for Information matrices_2
#########################################################

first_mix_num_hdf_11=function(i,m,n,nd,t)
{
  a=1
  first_mix_nf=first_mix_calc(t)[Cat_num_sim_hdf_11(i,a,m,n,nd)]
  return(first_mix_nf)
}

Prob_num_hdf_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_hdf_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hdf_11(i,a,m,n,nd),a,t)* 
                 (Sens_it[Cat_num_sim_hdf_11(i,a,m,n,nd),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_hdf_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hdf_11(i,a,m,n,nd),a,t)*
                 (1-Sens_it[Cat_num_sim_hdf_11(i,a,m,n,nd),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_hdf_11(i,a,m,n,nd),a,t)*(1-pr_diag(Cat_num_sim_hdf_11(i,a,m,n,nd),a,t)))
  }
  return(probab1)
}

Tprob_num_hdf_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_hdf_11(i,a,m,n,nd)),(Cat_num_sim_hdf_11(i,(a+1),m,n,nd))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_hdf_11(i,a,m,n,nd)),(Cat_num_sim_hdf_11(i,(a+1),m,n,nd))])
    
  }
  return(tprobab1)
}

primer_num_hdf_11=function(i,m,n,nd,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_hdf_11(i,m,n,nd,t)*Prob_num_hdf_11(i,a,m,n,nd,t)
  }
  else
  {
    a=1
    num=first_mix_num_hdf_11(i,m,n,nd,t)*Prob_num_hdf_11(i,a,m,n,nd,t)
    a=2
    repeat
    {
      num=num*Tprob_num_hdf_11(i,a-1,m,n,nd,t)*Prob_num_hdf_11(i,a,m,n,nd,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer2_num_hdf_11=function(i,n,nd,t)
{
  
  numh1=0
  m=1
  repeat
  {
    numh1=numh1+primer_num_hdf_11(i,m,n,nd,t)
    if(m==(dim(Hkdf11(i,n,nd))[[1]])){break}
    m=m+1
  }
  check2=numh1
  return(check2)
}

Joint_hdf_11=function(i,n,nd,t)
{
  j=(primer2_num_hdf_11(i,n,nd,t))/((primer2_den_h(i,t)))
  return(j)
}

#########################################################
# Similar set of computations for Information matrices_2
#########################################################

first_mix_num_hds_11=function(i,m,n,nd,t)
{
  a=1
  first_mix_ns=first_mix_calc(t)[Cat_num_sim_hds_11(i,a,m,n,nd)]
  return(first_mix_ns)
}

Prob_num_hds_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_hds_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hds_11(i,a,m,n,nd),a,t)* 
                 (Sens_it[Cat_num_sim_hds_11(i,a,m,n,nd),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_hds_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hds_11(i,a,m,n,nd),a,t)*
                 (1-Sens_it[Cat_num_sim_hds_11(i,a,m,n,nd),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_hds_11(i,a,m,n,nd),a,t)*(1-pr_diag(Cat_num_sim_hds_11(i,a,m,n,nd),a,t)))
  }
  return(probab1)
}

Tprob_num_hds_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_hds_11(i,a,m,n,nd)),(Cat_num_sim_hds_11(i,(a+1),m,n,nd))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_hds_11(i,a,m,n,nd)),(Cat_num_sim_hds_11(i,(a+1),m,n,nd))])
    
  }
  return(tprobab1)
}

primer_num_hds_11=function(i,m,n,nd,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_hds_11(i,m,n,nd,t)*Prob_num_hds_11(i,a,m,n,nd,t)
  }
  else
  {
    a=1
    num=first_mix_num_hds_11(i,m,n,nd,t)*Prob_num_hds_11(i,a,m,n,nd,t)
    a=2
    repeat
    {
      num=num*Tprob_num_hds_11(i,a-1,m,n,nd,t)*Prob_num_hds_11(i,a,m,n,nd,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer2_num_hds_11=function(i,n,nd,t)
{
  
  numh1=0
  m=1
  repeat
  {
    numh1=numh1+primer_num_hds_11(i,m,n,nd,t)
    if(m==(dim(Hkds11(i,n,nd))[[1]])){break}
    m=m+1
  }
  check2=numh1
  return(check2)
}

Joint_hds_11=function(i,n,nd,t)
{
  j=(primer2_num_hds_11(i,n,nd,t))/((primer2_den_h(i,t)))
  return(j)
}

#########################################################
# Similar set of computations for Information matrices_3
#########################################################

first_mix_num_hdw_11=function(i,m,n,nd,t)
{
  a=1
  first_mix_nw=first_mix_calc(t)[Cat_num_sim_hdw_11(i,a,m,n,nd)]
  return(first_mix_nw)
}

Prob_num_hdw_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==1)
    {
      probab1=((Diag[i,a])*  (C_Diag[i,a])*pr(i,Cat_num_sim_hdw_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hdw_11(i,a,m,n,nd),a,t)* 
                 (Sens_it[Cat_num_sim_hdw_11(i,a,m,n,nd),t]))
    }
    else
    {
      probab1=((Diag[i,a])*(1-C_Diag[i,a])*pr(i,Cat_num_sim_hdw_11(i,a,m,n,nd),a,t)*   pr_diag(Cat_num_sim_hdw_11(i,a,m,n,nd),a,t)*
                 (1-Sens_it[Cat_num_sim_hdw_11(i,a,m,n,nd),t]))
    }
  }
  else
  {
    probab1=((1-Diag[i,a])                *pr(i,Cat_num_sim_hdw_11(i,a,m,n,nd),a,t)*(1-pr_diag(Cat_num_sim_hdw_11(i,a,m,n,nd),a,t)))
  }
  return(probab1)
}

Tprob_num_hdw_11=function(i,a,m,n,nd,t)
{
  if(Diag[i,a]==1)
  {
    if(C_Diag[i,a]==0)
    {
      tprobab1= ((Diag[i,a])*(1-C_Diag[i,a])*TPM_calc(t)[(Cat_num_sim_hdw_11(i,a,m,n,nd)),(Cat_num_sim_hdw_11(i,(a+1),m,n,nd))])
    }
  }
  else
  {
    tprobab1=((1-Diag[i,a]) *TPM_calc(t)[(Cat_num_sim_hdw_11(i,a,m,n,nd)),(Cat_num_sim_hdw_11(i,(a+1),m,n,nd))])
    
  }
  return(tprobab1)
}

primer_num_hdw_11=function(i,m,n,nd,t)
{
  if(Stagescount(i)==1)
  {
    a=1
    num=first_mix_num_hdw_11(i,m,n,nd,t)*Prob_num_hdw_11(i,a,m,n,nd,t)
  }
  else
  {
    a=1
    num=first_mix_num_hdw_11(i,m,n,nd,t)*Prob_num_hdw_11(i,a,m,n,nd,t)
    a=2
    repeat
    {
      num=num*Tprob_num_hdw_11(i,a-1,m,n,nd,t)*Prob_num_hdw_11(i,a,m,n,nd,t)
      if(a==Stagescount(i)){break}
      a=a+1
    }
  }
  check=num
  return(check)
}

primer2_num_hdw_11=function(i,n,nd,t)
{
  
  numh1=0
  m=1
  repeat
  {
    numh1=numh1+primer_num_hdw_11(i,m,n,nd,t)
    if(m==(dim(Hkdw11(i,n,nd))[[1]])){break}
    m=m+1
  }
  check2=numh1
  return(check2)
}

Joint_hdw_11=function(i,n,nd,t)
{
  j=(primer2_num_hdw_11(i,n,nd,t))/((primer2_den_h(i,t)))
  return(j)
}

ptm <- proc.time()
CoVariance_2=function(i,a,ad,t,e,f)
{
  Co=0
  if(!is.na(Cat[i,a])&!is.na(Cat[i,ad]))
  {
  if((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,ad]==3|Cat[i,ad]==4))
  {
    if((e<3)&(f<3))
    {
      Co=((-1)^((e+f)%%2))*(Joint_hd_11(i,a,ad,t)-Expzee(i,a,t)*Expzee(i,ad,t))
      
    }
    else
    {
      Co=0
    }
  }
  else
  {
    Co=0
  }
  }
  return(Co)
}
proc.time() - ptm

CoVariance_3=function(i,a,ad,t,e,f,g)
{
  CoVari=0
  if(!is.na(Cat[i,a])&!is.na(Cat[i,ad])&!is.na(Cat[i,ad+1]))
  {
  if((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,ad]==3|Cat[i,ad]==4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4))
  {
   if(e==1&f==1&g==1)
   {
    CoVari=Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)
   }
   if(e==1&f==1&g==2)
   {
    Covari=Joint_hd_11(i,a,ad,t)-Expzee(i,a,t)*Expzee(i,ad,t)+Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdf_11(i,a,ad,t)
   }
   if(e==1&f==2&g==1)
   {
     Covari=Joint_hd_11(i,a,ad+1,t)-Expzee(i,a,t)*Expzee(i,ad+1,t)+Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdf_11(i,a,ad,t)
   }
   if(e==2&f==1&g==1)
   {
     CoVari=Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdf_11(i,a,ad,t)
   }
   if(e==1&f==2&g==2)
   {
     CoVari=Expzee(i,a,t)*Expzee(i,ad,t)-Joint_hd_11(i,a,ad,t)+Expzee(i,a,t)*Expzee(i,ad+1,t)-Joint_hd_11(i,a,ad+1,t)+Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)
   }
   if(e==2&f==1&g==2)
   {
    CoVari=Joint_hd_11(i,a,ad,t)-Expzee(i,a,t)*Expzee(i,ad,t)+Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)
   }
   if(e==2&f==2&g==1)
   {
     CoVari=Joint_hd_11(i,a,ad+1,t)-Expzee(i,a,t)*Expzee(i,ad+1,t)+Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)
   }
   if(e==2&f==2&g==2)
   {
     CoVari=Joint_hd_11(i,a,ad+1,t)-Expzee(i,a,t)*Expzee(i,ad+1,t)+Joint_hd_11(i,a,ad,t)-Expzee(i,a,t)*Expzee(i,ad,t)-Joint_hdf_11(i,a,ad,t)+Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)
   }
  }
  if((Cat[i,ad]!=3&Cat[i,ad]!=4)&(Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4))
  {
    if(Cat[i,ad]==f)
    {
      CoVari=CoVariance_2(i,a,ad+1,t,e,g)
    }
  }
  if((Cat[i,ad+1]!=3&Cat[i,ad+1]!=4)&(Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,ad]==3|Cat[i,ad]==4))
  {
    if(Cat[i,ad+1]==g)
    {
      CoVari=CoVariance_2(i,a,ad,t,e,f)
    }
  }
  }
  
  return(CoVari)
    
}
#CoVariance_c4(i,k,kd,iter,e,f,g,h)
# a=k
# ad=kd

CoVariance_4=function(i,a,ad,t,e,f,g,h)
{
  CoVari=0
  if(!is.na(Cat[i,a])&!is.na(Cat[i,ad])&!is.na(Cat[i,ad+1])&!is.na(Cat[i,a+1]))
  {
  if((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,a+1]==3|Cat[i,a+1]==4)&(Cat[i,ad]==3|Cat[i,ad]==4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4))
  {
    if(e==1&f==1&g==1&h==1)
    {
      CoVari=Joint_hdw_11(i,a,ad,t)-(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==1&f==1&g==1&h==2)
    {
     CoVari=Joint_hdf_11(i,ad,a,t)-Expzee(i,ad,t)*Joint_hd_11(i,a,a+1,t)-Joint_hdw_11(i,a,ad,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==1&f==1&g==2&h==1)
    {
      CoVari=Joint_hdf_11(i,ad+1,a,t)-Expzee(i,ad+1,t)*Joint_hd_11(i,a,a+1,t)-Joint_hdw_11(i,a,ad,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==1&f==2&g==1&h==1)
    {
      CoVari=Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdw_11(i,a,ad,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==2&f==1&g==1&h==1)
    {
      CoVari=Joint_hdf_11(i,a+1,ad,t)-Expzee(i,a+1,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdw_11(i,a,ad,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==1&f==1&g==2&h==2)
    {
      CoVari=Expzee(i,ad,t)*Joint_hd_11(i,a,a+1,t)-Joint_hdf_11(i,ad,a,t)+Expzee(i,ad+1,t)*Joint_hd_11(i,a,a+1,t)-Joint_hdf_11(i,ad+1,a,t)+Joint_hdw_11(i,a,ad,t)-(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==2&f==2&g==1&h==1)
    {
      CoVari=Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdf_11(i,a,ad,t)+Expzee(i,a+1,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdf_11(i,a+1,ad,t)+Joint_hdw_11(i,a,ad,t)-(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==2&f==1&g==1&h==2)
    {
      CoVari=Joint_hd_11(i,a+1,ad,t)-Expzee(i,a+1,t)*Expzee(i,ad,t)+Joint_hd_11(i,a,a+1,t)*Expzee(i,ad,t)-Joint_hdf_11(i,ad,a,t)+Joint_hd_11(i,ad,ad+1,t)*Expzee(i,a+1,t)-Joint_hdf_11(i,a+1,ad,t)+Joint_hdw_11(i,a,ad,t)-(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==1&f==2&g==2&h==1)
    {
      CoVari=Joint_hd_11(i,a,ad+1,t)-Expzee(i,a,t)*Expzee(i,ad+1,t)+Joint_hd_11(i,a,a+1,t)*Expzee(i,ad+1,t)-Joint_hdf_11(i,ad+1,a,t)+Joint_hd_11(i,ad,ad+1,t)*Expzee(i,a,t)-Joint_hdf_11(i,a,ad,t)+Joint_hdw_11(i,a,ad,t)-(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==2&f==1&g==2&h==1)
    {
      CoVari=Joint_hd_11(i,a+1,ad+1,t)-Expzee(i,a+1,t)*Expzee(i,ad+1,t)+Joint_hd_11(i,a,a+1,t)*Expzee(i,ad+1,t)-Joint_hdf_11(i,ad+1,a,t)+Joint_hd_11(i,ad,ad+1,t)*Expzee(i,a+1,t)-Joint_hdf_11(i,a+1,ad,t)+Joint_hdw_11(i,a,ad,t)-(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==1&f==2&g==1&h==2)
    {
      CoVari=Joint_hd_11(i,a,ad,t)-Expzee(i,a,t)*Expzee(i,ad,t)+Joint_hd_11(i,a,a+1,t)*Expzee(i,ad,t)-Joint_hdf_11(i,ad,a,t)+Joint_hd_11(i,ad,ad+1,t)*Expzee(i,a,t)-Joint_hdf_11(i,a,ad,t)+Joint_hdw_11(i,a,ad,t)-(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))
    }
    if(e==2&f==2&g==2&h==1)
    {
      CoVari=Expzee(i,a,t)*Expzee(i,ad+1,t)-Joint_hd_11(i,a,ad+1,t)+Expzee(i,a+1,t)*Expzee(i,ad+1,t)-Joint_hd_11(i,a+1,ad+1,t)+Joint_hdf_11(i,ad+1,a,t)-Expzee(i,ad+1,t)*Joint_hd_11(i,a,a+1,t)+Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)+Joint_hdf_11(i,a+1,ad,t)-Expzee(i,a+1,t)*Joint_hd_11(i,ad,ad+1,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))-Joint_hdw_11(i,a,ad,t)
    }
    
    if(e==2&f==2&g==1&h==2)
    {
      CoVari=Expzee(i,a,t)*Expzee(i,ad,t)-Joint_hd_11(i,a,ad,t)+Expzee(i,a+1,t)*Expzee(i,ad,t)-Joint_hd_11(i,a+1,ad,t)+Joint_hdf_11(i,ad,a,t)-Expzee(i,ad,t)*Joint_hd_11(i,a,a+1,t)+Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)+Joint_hdf_11(i,a+1,ad,t)-Expzee(i,a+1,t)*Joint_hd_11(i,ad,ad+1,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))-Joint_hdw_11(i,a,ad,t)
    }
    if(e==2&f==1&g==2&h==2)
    {
      CoVari=Expzee(i,ad,t)*Expzee(i,a+1,t)-Joint_hd_11(i,ad,a+1,t)+Expzee(i,ad+1,t)*Expzee(i,a+1,t)-Joint_hd_11(i,ad+1,a+1,t)+Joint_hdf_11(i,a+1,ad,t)-Expzee(i,a+1,t)*Joint_hd_11(i,ad,ad+1,t)+Joint_hdf_11(i,ad,a,t)-Expzee(i,ad,t)*Joint_hd_11(i,a,a+1,t)+Joint_hdf_11(i,ad+1,a,t)-Expzee(i,ad+1,t)*Joint_hd_11(i,a,a+1,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))-Joint_hdw_11(i,a,ad,t)
    }
    if(e==1&f==2&g==2&h==2)
    {
      CoVari=Expzee(i,ad,t)*Expzee(i,a,t)-Joint_hd_11(i,ad,a,t)+Expzee(i,ad+1,t)*Expzee(i,a,t)-Joint_hd_11(i,ad+1,a,t)+Joint_hdf_11(i,a,ad,t)-Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)+Joint_hdf_11(i,ad,a,t)-Expzee(i,ad,t)*Joint_hd_11(i,a,a+1,t)+Joint_hdf_11(i,ad+1,a,t)-Expzee(i,ad+1,t)*Joint_hd_11(i,a,a+1,t)+(Joint_hd_11(i,a,(a+1),t))*(Joint_hd_11(i,ad,(ad+1),t))-Joint_hdw_11(i,a,ad,t)
    }
    if(e==2&f==2&g==2&h==2)
    {
      CoVari=Joint_hd_11(i,a,ad,t)-Expzee(i,a,t)*Expzee(i,ad,t)+Joint_hd_11(i,a,ad+1,t)-Expzee(i,a,t)*Expzee(i,ad+1,t)+Joint_hd_11(i,a+1,ad,t)-Expzee(i,a+1,t)*Expzee(i,ad,t)+Joint_hd_11(i,a+1,ad+1,t)-Expzee(i,a+1,t)*Expzee(i,ad+1,t)+Expzee(i,a,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdf_11(i,a,ad,t)+Expzee(i,a+1,t)*Joint_hd_11(i,ad,ad+1,t)-Joint_hdf_11(i,a+1,ad,t)+Expzee(i,ad,t)*Joint_hd_11(i,a,a+1,t)-Joint_hdf_11(i,ad,a,t)+Expzee(i,ad+1,t)*Joint_hd_11(i,a,a+1,t)-Joint_hdf_11(i,ad+1,a,t)+Joint_hdw_11(i,a,ad,t)-Joint_hd_11(i,a,a+1,t)*Joint_hd_11(i,ad,ad+1,t)
    }
  }
  if(((Cat[i,a]!=3&Cat[i,a]!=4)&(Cat[i,a+1]==3|Cat[i,a+1]==4)&(Cat[i,ad]==3|Cat[i,ad]==4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4)))
  {
    if(Cat[i,a]==e)
    {
      CoVari=CoVariance_3(i,a+1,ad,t,f,g,h)
    }
  }
  if(((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,a+1]!=3&Cat[i,a+1]!=4)&(Cat[i,ad]==3|Cat[i,ad]==4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4)))
  {
    if(Cat[i,a+1]==f)
    {
      CoVari=CoVariance_3(i,a,ad,t,e,g,h)
    }
  }
  if(((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,a+1]==3|Cat[i,a+1]==4)&(Cat[i,ad]!=3&Cat[i,ad]!=4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4)))
  {
    if(Cat[i,ad]==g)
    {
      CoVari=CoVariance_3(i,ad+1,a,t,h,e,f)
    }
  }
  if(((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,a+1]==3|Cat[i,a+1]==4)&(Cat[i,ad]==3|Cat[i,ad]==4)&(Cat[i,ad+1]!=3&Cat[i,ad+1]!=4)))
  {
    if(Cat[i,ad+1]==h)
    {
      CoVari=CoVariance_3(i,ad,a,t,g,e,f)
    }
  }
  if(((Cat[i,a]!=3&Cat[i,a]!=4)&(Cat[i,a+1]==3|Cat[i,a+1]==4)&(Cat[i,ad]!=3&Cat[i,ad]!=4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4)))
  {
    if(Cat[i,a]==e&Cat[i,ad]==g)
    {
      CoVari=CoVariance_2(i,a+1,ad+1,t,f,h)
    }
  }
  
  if(((Cat[i,a]!=3&Cat[i,a]!=4)&(Cat[i,a+1]==3|Cat[i,a+1]==4)&(Cat[i,ad]==3|Cat[i,ad]==4)&(Cat[i,ad+1]!=3&Cat[i,ad+1]!=4)))
  {
    if(Cat[i,a]==e&Cat[i,ad+1]==h)
    {
      CoVari=CoVariance_2(i,a+1,ad,t,f,g)
    }
  }
  if(((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,a+1]!=3&Cat[i,a+1]!=4)&(Cat[i,ad]!=3&Cat[i,ad]!=4)&(Cat[i,ad+1]==3|Cat[i,ad+1]==4)))
  {
    if(Cat[i,a+1]==f&Cat[i,ad]==g)
    {
      CoVari=CoVariance_2(i,a,ad+1,t,e,h)
    }
  }
  if(((Cat[i,a]==3|Cat[i,a]==4)&(Cat[i,a+1]!=3&Cat[i,a+1]!=4)&(Cat[i,ad]==3|Cat[i,ad]==4)&(Cat[i,ad+1]!=3&Cat[i,ad+1]!=4)))
  {
    if(Cat[i,a+1]==f&Cat[i,ad+1]==h)
    {
      CoVari=CoVariance_2(i,a,ad,t,e,g)
    }
  }
  }
    
return(CoVari)
}


proc.time() - ptm