##############################
# Information Matrices
##############################

##Complete Matrices


Constt0=Category-1
Constt1=Constt0+Category
Constt2=Constt1+Category
Constt3=Constt2+Category
Constt4=Constt3+Category*(Category-1)
Number_Parameters_Estimated=Constt4

I_C=array(0, dim=c(Number_Parameters_Estimated,Number_Parameters_Estimated))

# Computation of I_C

k=1
I_C[1:Constt0, 1:Constt0]=(sum(Observed(Numb[,Category,k])))/((first_mix[Category])^2)

for(e in 1:(Category-1))
{
  k=1
  I_C[e,e]= (sum(Observed(Numb[,e,k]))*((1/first_mix[e])^2))+
    (sum(Observed(Numb[,Category,k]))*((1/first_mix[Category])^2))
}

p=0
for(e in 1:Category)
{
    j=Constt0+e
    for(k in 1:Stages_Interval)
      {
      p=p+sum(Observed(Numb[,e,k]*Diag[,k]))
    }
    I_C[j,j]=p/((mu_d[e])^2)
  }
q=0

for(e in 1:Category)
{
  j=Constt1+e
  for(k in 1:Stages_Interval)      
  {
    q=q+sum(Observed((Numb[,e,k])*(1-Diag[,k])))
    
  }
  I_C[j,j]=q/((mu_s[e])^2)
}
r=0
m=0
for(e in 1:Category)
{
  j=Constt2+e
  for(k in 1:Stages_Interval)      
  {
    r=r+sum(Observed((Numb[,e,k])*(Diag[,k])*(C_Diag[,k])))
    m=m+sum(Observed((Numb[,e,k])*(Diag[,k])*(1-C_Diag[,k])))
    
}
  I_C[j,j]=r/((Sens[e])^2)+m/(1-Sens[e])^2
}


for(e in 1:Category)
{
  for(f in 1:(Category-1))
  {
    j=Constt3+((e-1)*(Category-1))+f
    I_C[j,j]=(sum(Observed(TPM_numb[e,f])))/((TPM[e,f])^2)+
      (sum(Observed(TPM_numb[e,Category])))/((TPM[e,Category])^2)
  }
}
for(e in 1:Category)
{
for(f in 1:(Category-1))
{
  for(g in 1:(Category-1))
  {
  q=Constt3+(e-1)*(Category-1)+f
  r=Constt3+(e-1)*(Category-1)+g
  if(q!=r)
  {
  I_C[q,r]=(sum(Observed(TPM_numb[e,Category])))/((TPM[e,Category])^2)
  }
}
}
}




write.csv(I_C,file="parameters/I_C.csv")
