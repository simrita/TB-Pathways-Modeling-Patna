
W_Diag=Diag[,1:Stages_Interval]-C_Diag[,1:Stages_Interval]
No_Diag=1-Diag[,1:Stages_Interval]
wrong_waits<-NULL
switch_waits<-NULL
correct_waits<-NULL

wrong_waits=vector()
switch_waits=vector()
correct_waits=vector()
for(i in 1:n_obs){
  for(k in 1:(Stages_Interval-1)){
    if(!is.na(C_Diag[i,k])){
      if(C_Diag[i,k]==1){
        correct_waits<-append(correct_waits,Waiting_time[i,k])
      }
    }
  }
}

for(i in 1:n_obs){
  for(k in 1:(Stages_Interval-1)){
    if(!is.na(Cat[i,k])){
      if(W_Diag[i,k]==1){
        wrong_waits<-append(wrong_waits,Waiting_time[i,k])
      }
    }
  }
}
for(i in 1:n_obs){
  for(k in 1:(Stages_Interval-1)){
    if(!is.na(Cat[i,k])){
      if(No_Diag[i,k]==1){
        switch_waits<-append(switch_waits,Waiting_time[i,k])
      }
    }
  }
}

waits<-c(wrong_waits,switch_waits) ########## for sampling wait times
waits<-Observed(waits)



library(MASS)
var_matrix_base=Inverse
mean_base=c(first_mix[1:3],mu_d[1:2],mu_d[4],mu_s,Sens[2],TPM[1,1:2],TPM[2,2],TPM[3,1:3])
mean_base<-as.numeric(mean_base)

generate_parameters=function(mean_vector,var_matrix){
  temp=mvrnorm(n=1,mu=mean_vector,Sigma=var_matrix)
  for(i in 1:length(temp)){
    
    if(temp[i]<0){
      temp[i]<-0
    }
  }
  for(i in 12:17){
    if(temp[i]>1){
      temp[i]=1
    }
  }

  


  
  parameter_vector=vector()
  parameter_vector[1:3]=temp[1:3]
  parameter_vector[4]=max(1-sum(temp[1:3]),0)
  parameter_vector[5:6]=temp[4:5]
  parameter_vector[7]=0
  parameter_vector[8]=temp[6]
  parameter_vector[9:12]=temp[7:10]
  parameter_vector[13]=1
  parameter_vector[14]=temp[11]
  parameter_vector[15]=0
  parameter_vector[16]=1
  parameter_vector[17:18]=temp[12:13]
  parameter_vector[19]=0
  parameter_vector[20]=max(1-sum(temp[12:13]),0)
  parameter_vector[21]=0
  parameter_vector[22]=temp[14]
  parameter_vector[23]=0
  parameter_vector[24]=max(1-temp[14],0)
  parameter_vector[25:27]=temp[15:17]
  
  parameter_vector[28]=max(1-sum(temp[15:17]),0)
  parameter_vector[29:31]=0
 
  parameter_vector[32]=1
  return(parameter_vector)
}
###### distribution of the diagnostic delay and number of stages of the original data
diagnostic_delay_original=vector()
delay_original=data.frame(dim=c(n_obs,Stages_Interval))
waits_delay_original=data.frame(dim=c(n_obs,Stages_Interval))
delta_diff_wrong_original=data.frame(dim=c(n_obs,Stages_Interval))
for(i in 1:n_obs){
  for(k in 1:Stages_Interval){
    if(!is.na(Diag[i,k])){
      if(C_Diag[i,k]==1){
        delay_original[i,k]=Dur[i,k]
      }
      if(Diag[i,k]==0){
        delay_original[i,k]=Dur[i,k]+Waiting_time[i,k]
      }
      if(Diag[i,k]==1 & C_Diag[i,k]==0){
        delay_original[i,k]=Dur[i,k]+Waiting_time[i,k]+Delta_diffe_wrong[i,k]
      }
    }
  }
}
for(i in 1:n_obs){
  diagnostic_delay_original[i]=sum(Observed(delay_original[i,]))
}
pathway_stages_original=vector()
for(i in 1:n_obs){
  pathway_stages_original[i]=Stagescount(i)
}
diagnostic_delay_sim=vector()
pathway_stages_sim=vector()
simulations=10000
stages=5
cat=array(dim=c(simulations,stages))### indicates the category of each patient at each stage
diag=array(dim=c(simulations,stages))### indidicates whether patient was diagnosed or not
cdiag=array(dim=c(simulations,stages))#### indicates whether patient was correctly diagnosed or not
wdiag=diag-cdiag ##### indicates whether patient was wrongly diagnosed or not
swtch<- 1-diag ##### indicates whether patient switched before diagnosis or not
d_dur<-array(dim=c(simulations,stages))
s_dur<-array(dim=c(simulations,stages))
dur<-array(dim=c(simulations,stages))

   for(j in 1:simulations){
    para=generate_parameters(mean_base,var_matrix_base)
    sim_first_mix=para[1:4]
    sim_mu_d=para[5:8]
    sim_mu_s=para[9:12]
    sim_sens=para[13:16]
    sim_TPM=array(dim=c(Category,Category))
    sim_TPM[1,]=para[17:20]
    sim_TPM[2,]=para[21:24]
    sim_TPM[3,]=para[25:28]
    sim_TPM[4,]=para[29:32]
    
      #patient was diagnosed at the first stage or not
      a=runif(1)
      b=runif(1)
      c=runif(1)
      d=runif(1)
      e=runif(1)
      f=runif(1)
      g=runif(1)
      h=runif(1)
      m=runif(1)
      n=runif(1)
      #### to allocate the category for the first stage
      for(k in 1:4){
        if (a <= cumsum(sim_first_mix)[k])
        {
          cat[j,1]<-k
          break
        }
        
      }
      if(cat[j,1]!=3){
        d_dur[j,1]<-rexp(1,rate=sim_mu_d[cat[j,1]]+0.0000000000000000000000000000000000001)
      }else{
        d_dur[j,1]<-9999999999999999999999999999999999
      }
      s_dur[j,1]<-rexp(1,rate=sim_mu_s[cat[j,1]])
      dur[j,1]<-min(d_dur[j,1],s_dur[j,1])
      if(d_dur[j,1] < s_dur[j,1])
      {
        diag[j,1]<-1
      } else{
        diag[j,1]<-0
      }
      ##### to decide whether diagnosis was correct or not
      if(diag[j,1]==1 && b <= sim_sens[cat[j,1]])
      {
        cdiag[j,1]<-1
      }else{
        cdiag[j,1]<-0
      }
      wdiag[j,1]<-diag[j,1]-cdiag[j,1]
      swtch[j,1]<-1-diag[j,1]
      #### To find the next category i.e cat[j,2]
      ####if(cdiag[j,1] == 1){break}
      
      if(cdiag[j,1]==0){
        
        for(k in 1:4){
          if (c <= cumsum(sim_TPM[cat[j,1],])[k])
          {
            cat[j,2]<-k
            break
          }
        }
        
        
        
        ### to decide what happened at the second stage
        
        
        if(cat[j,2]!=3){
          d_dur[j,2]<-rexp(1,rate=sim_mu_d[cat[j,2]])
        }else{
          d_dur[j,2]=999999999999999999999999999999999
        }
        s_dur[j,2]<-rexp(1,rate=sim_mu_s[cat[j,2]])
        dur[j,2]<-min(d_dur[j,2],s_dur[j,2])
        if(d_dur[j,2] < s_dur[j,2])
        {
          diag[j,2]<-1
        } else{
          diag[j,2]<-0
        }
        ##### to decide whether diagnosis was correct or not
        if(diag[j,2]==1 && d <= sim_sens[cat[j,2]])
        {
          cdiag[j,2]<-1
        }else{
          cdiag[j,2]<-0
        }
        wdiag[j,2]<-diag[j,2]-cdiag[j,2]
        swtch[j,2]<-1-diag[j,2]
        #### To find the next category i.e cat[j,3]
        
        if(cdiag[j,2]==0){
          
          for(k in 1:4){
            if (e <= cumsum(sim_TPM[cat[j,2],])[k])
            {
              cat[j,3]<-k
              break
              
            }
            
          }
          
          
          ######## to decide what happened at 3rd stage
          
          if(cat[j,3]!=3){
            d_dur[j,3]<-rexp(1,rate=sim_mu_d[cat[j,3]])
          }else{
            d_dur[j,3]<-99999999999999999999999999999999
          }
          s_dur[j,3]<-rexp(1,rate=sim_mu_s[cat[j,3]])
          dur[j,3]<-min(d_dur[j,3],s_dur[j,3])
          if(d_dur[j,3] < s_dur[j,3])
          {
            diag[j,3]<-1
          } else{
            diag[j,3]<-0
          }
          ##### to decide whether diagnosis was correct or not
          if(diag[j,3]==1 && f <= sim_sens[cat[j,3]])
          {
            cdiag[j,3]<-1
          }else{
            cdiag[j,3]<-0
          }
          wdiag[j,3]<-diag[j,3]-cdiag[j,3]
          swtch[j,3]<-1-diag[j,3]
          #### To find the next category i.e cat[j,4]
          
          if(cdiag[j,3]==0){
            
            for(k in 1:4){
              if (g <= cumsum(sim_TPM[cat[j,3],])[k])
              {
                cat[j,4]<-k
                break
                
              }
              
            }
            
            
            ######## to decide what happened at 4th stage
            
            if(cat[j,4]!=3){
              d_dur[j,4]<-rexp(1,rate=sim_mu_d[cat[j,4]])
            }else{
              d_dur[j,4]<-999999999999999999999999999999
            }
            s_dur[j,4]<-rexp(1,rate=sim_mu_s[cat[j,4]])
            dur[j,4]<-min(d_dur[j,4],s_dur[j,4])
            if(d_dur[j,4] < s_dur[j,4])
            {
              diag[j,4]<-1
            } else{
              diag[j,4]<-0
            }
            ##### to decide whether diagnosis was correct or not
            if(diag[j,4]==1 && h <= sim_sens[cat[j,4]])
            {
              cdiag[j,4]<-1
            }else{
              cdiag[j,4]<-0
            }
            wdiag[j,4]<-diag[j,4]-cdiag[j,4]
            swtch[j,4]<-1-diag[j,4]
            
            if(cdiag[j,4]==0){
              
              for(k in 1:4){
                if (m <= cumsum(sim_TPM[cat[j,4],])[k])
                {
                  cat[j,5]<-k
                  break
                  
                }
                
              }
              
              if(cat[j,5]!=3)
              {
                d_dur[j,5]<-rexp(1,rate=sim_mu_d[cat[j,5]])
              }else{
                d_dur[j,5]<-99999999999999999999999999999
              }
              s_dur[j,5]<-rexp(1,rate=sim_mu_s[cat[j,5]])
              dur[j,5]<-min(d_dur[j,5],s_dur[j,5])
              if(d_dur[j,5] < s_dur[j,5])
              {
                diag[j,5]<-1
              } else{
                diag[j,5]<-0
              }
              ##### to decide whether diagnosis was correct or not
              if(diag[j,5]==1 && n <= sim_sens[cat[j,5]])
              {
                cdiag[j,5]<-1
              }else{
                cdiag[j,5]<-0
              }
              wdiag[j,5]<-diag[j,5]-cdiag[j,5]
              swtch[j,5]<-1-diag[j,5]
              
              
            }
          }
        }
      }
      
    }
    
    
    
    diagnostic_delay_sim=vector()
    pathway_stages_sim=vector()
    
    delay_temp=data.frame(dim=c(simulations,stages))

    
    for(p in 1:simulations){
      for(q in 1:stages){
        if(!is.na(cat[p,q])){
          if(cdiag[p,q]==1){
            delay_temp[p,q]=dur[p,q]
          }
          if(diag[p,q]==0){
            delay_temp[p,q]=dur[p,q]+sample(waits,1)
          }
          if(diag[p,q]==1 & cdiag[p,q]==0){
            delay_temp[p,q]=dur[p,q]+sample(waits,1)+sample(Delta_diff_wrong,1)
            
          }
        }
      }
    }
    for(j in 1:simulations){
      diagnostic_delay_sim[j]=sum(Observed(delay_temp[j,]))
    }
    
    for(j in 1:simulations){
      pathway_stages_sim[j]=length(Observed(cat[j,]))
    }
    
    