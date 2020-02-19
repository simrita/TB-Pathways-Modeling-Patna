# This fragment of code computes the intervals and draws the plots of the estiamtes with their error intervals.

Order=array(dim=Number_Parameters_Estimated)
z=array(dim=Number_Parameters_Estimated)

 I=I_C-I_M

for(i in 1:length(z))
{
  Order[i]=i
  if((diag(I)[i]==0)|(diag(I)[i]<0)|(diag(I)[i]==Inf)|(is.na(diag(I))[i]=="TRUE"))
  {
    z[i]=1
  }
  else
    
  {
    z[i]=0
  }
}

r=Order[z==1]
I_1=I[-r,-r]

write.csv(I_1,file="parameters/I.csv")

Inverse=qr.solve(I_1)
SE=diag(qr.solve(I_1)^0.5)

error=array(dim=Constt4)
j=1
for(i in 1:Number_Parameters_Estimated)
{
  if(z[i]==1)
  {
    error[i]=0
  }
  else
  {
    error[i]=SE[j]
    j=j+1
  }
}


Const0=Category
Const1=Const0+(Category)
Const2=Const1+(Category)
Const3=Const2+(Category)
Const4=Const3+(Category*Category)

errorf=array(dim=Const4)

errorf[1:Constt0]=error[1:Constt0]
errorf[Const0]=sqrt(sum((error[1:Constt0])^2))

errorf[(Const0+1):Const1]=error[(Constt0+1):Constt1]
errorf[(Const1+1):Const2]=error[(Constt1+1):Constt2]
errorf[(Const2+1):Const3]=error[(Constt2+1):Constt3]

errorf[Const3+1:3]=error[Constt3+1:3]
errorf[Const3+4]=sqrt(sum((errorf[Const3+1:3])^2))

p=Const3+4
pd=Constt3+3

errorf[p+1:3]=error[pd+1:3]
errorf[p+4]=sqrt(sum((errorf[p+1:3])^2))

p=p+4
pd=pd+3

errorf[p+1:3]=error[pd+1:3]
errorf[p+4]=sqrt(sum((errorf[p+1:3])^2))

p=p+4
pd=pd+3

errorf[p+1:3]=error[pd+1:3]
errorf[p+4]=sqrt(sum((errorf[p+1:3])^2))

write.csv(errorf,"errors/error.csv")

first_mix_e=array(dim=Category)
mu_d_e=vector()
mu_s_e=vector()
Sens_e=vector()
TPM_e=array(dim=c(Category, Category))


for(e in 1:Category)
{
  first_mix_e[e]=errorf[e]
}

write.csv(first_mix_e,"errors/first_mix_e.csv")

p=Const0

for(e in 1:Category)
{
mu_d_e[e]=errorf[p+e]
}

write.csv( mu_d_e,"errors/mu_d_e.csv")

p=Const1

for(e in 1:Category)
{
  
  mu_s_e[e]=errorf[p+e]
  
}

write.csv( mu_s_e, "errors/mu_s_e.csv")

p=Const2

for(e in 1:Category)
{
  
Sens_e[e]=errorf[p+e]
  
}

write.csv( Sens_e, "errors/Sens.csv")

p=Const3

for(e in 1:Category)
{
  for(f in 1:Category)
  {
    z=((e-1)*Category)+f
    TPM_e[e,f]=errorf[p+z]
  }
}

write.csv( TPM_e, "errors/TPM_e.csv")


########################################################
#COnfidence Intervals for given error 
########################################################

CI=function(x,error)
{
  zvalue1=1.96
  zvalue2=1.65
  upperlimit=x+(zvalue1*error)
  lowerlimit=x-(zvalue1*error)
  b2=x-(zvalue1*error)
  b1=x-(zvalue2*error)
  y=x/2
  
 if(!is.na(y))
{
if(y<=b2)
{
lowerlimit=b2
}
if(y>b2 & y<=b1)
{
lowerlimit=y
}
if(y>b1)
{
lowerlimit=max(b1,0)
}
}
  
 return(list(lowerlimit, upperlimit))
 }

######################
# Intervals
######################

first_mix_l=array(dim=c(Category,2))
mu_d_l=array(dim=c(Category,2))
mu_s_l=array(dim=c(Category,2))
Sens_l=array(dim=c(Category,2))
TPM_l=array(dim=c(Category, Category,2))


for(e in 1:Category)
{
  first_mix_l[e,1]=CI(first_mix[e],first_mix_e[e])[[1]]
  first_mix_l[e,2]=CI(first_mix[e],first_mix_e[e])[[2]]
}
first_mix_complete=data.frame(first_mix,first_mix_l)

write.csv(first_mix_complete, "confidence_intervals/first_mix_complete.csv")



for(e in 1:Category)
{

  
    mu_d_l[e,1]=CI(mu_d[e], mu_d_e[e])[[1]]
    mu_d_l[e,2]=CI(mu_d[e], mu_d_e[e])[[2]]
  
}

mu_d_complete=data.frame(mu_d,mu_d_l)
write.csv(mu_d_complete,"confidence_intervals/mu_d_complete.csv")
for(e in 1:Category)
{
  
  
  mu_s_l[e,1]=CI(mu_s[e], mu_s_e[e])[[1]]
  mu_s_l[e,2]=CI(mu_s[e], mu_s_e[e])[[2]]
  
}

mu_s_complete=data.frame(mu_s,mu_s_l)
write.csv(mu_s_complete,"confidence_intervals/mu_s_complete.csv")


for(e in 1:Category)
{
  
  
  Sens_l[e,1]=CI(Sens[e], Sens_e[e])[[1]]
  Sens_l[e,2]=CI(Sens[e], Sens_e[e])[[2]]
  
}

Sens_complete=data.frame(Sens,Sens_l)
write.csv(Sens_complete,"confidence_intervals/Sens_complete.csv")


for(e in 1:Category)
{
  for(f in 1:Category)
  {
    TPM_l[e,f,1]=CI(TPM[e,f], TPM_e[e,f])[[1]]
    TPM_l[e,f,2]=CI(TPM[e,f], TPM_e[e,f])[[2]]
  }
}

TPM_complete=data.frame(TPM,TPM_l)
write.csv(TPM_complete,"confidence_intervals/TPM_complete.csv")

