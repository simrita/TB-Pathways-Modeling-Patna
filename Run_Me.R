#### Before running this script ensure that all the excel files  and R scripts are placed in your R woking directory which can be checked by getwd()
#### The script creates the folders : parameters , errors and confidence_intervals in the R working directory.This folder contains the estimated parameters, standard errors and confidence intervals
# #### The script may take a couple of hours to run.
# 
# dir.create("parameters")
# dir.create("errors")
# dir.create("confidence_intervals")
source("Step-1-Data-Entry.r")
source("Step-2-Data-Extraction.r")
source("Step-3-Data-Extraction.r")
source("Step-4-Defining-Functions.r")
source("Step-5-EM.r")
write.csv(first_mix_it[,t+1], "parameters/First_mix.csv")
write.csv(1/mu_d_it[,t+1], "parameters/T_d.csv")
write.csv(1/mu_s_it[,t+1], "parameters/T_s.csv")
write.csv(Sens_it[,t+1], "parameters/Sensitivity.csv")

write.csv(TPM_it[,,t+1], "parameters/TPM.csv")


write.csv(Num[,1:Category,1,t], "parameters/Num_1.csv")
write.csv(Num[,1:Category,2,t], "parameters/Num_2.csv")
write.csv(Num[,1:Category,3,t], "parameters/Num_3.csv")
write.csv(Num[,1:Category,4,t], "parameters/Num_4.csv")
write.csv(Num[,1:Category,5,t], "parameters/Num_5.csv")

write.csv(W_TPM_num[,,t+1], "parameters/W_TPM_num.csv")
write.csv(S_TPM_num[,,t+1], "parameters/S_TPM_num.csv")

write.csv(Expz[,,t],"parameters/Expz.csv")

write.csv(TPM_num[,,t+1],"parameters/TPM_num.csv")
source("Step6-Reading parameters.r")
source("Step7-Ic.r")
source("Step8-Im.r")
source("Step 9 - Confidence intervals.r")

proc.time() - ptm
