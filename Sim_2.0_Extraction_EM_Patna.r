#############################################################
#This fragment of code extracts the observables from care pathway data.
#############################################################

#############################################################
#Additional Assumptions
###############################
#Initially there are five Category of Providers: LTFQs, FQs, Chemists, Public Facilities and Private Providers 
#with Unknown Qualifications

#Providers with MBBS as qualification are classified as Private FQs

#Providers with unknown qualifications have to be reclassified into LTFQs or FQs

#Reason to switch from the provider classiffied as pre/post diagnosis and provider-driven and patient driven

#Post diagnosis, three estimates are computed:  Delta-same(correct), Delta-different(correct/wrong) and 
#Delta-confirmation

#In the excel sheet replace entries as "medicines provided by health post don't work" with 
#"medicines provided by health post don't work"

ptm <- proc.time()

Stages_Estimation_1.0_TB = Stages_Estimation_2.0_TB= 4
# Stages_Estimation_2.0_TB = 7
Number_Category_Providers=5

#################################################################
#(Truncating patients with more than certain stages of consultation) 
#################################################################

#TB_Pat_Pathway=subset(TB_Pat_Pathway_6.0, TB_Pat_Pathway_6.0$D4_A4=="#NULL!")
# TB_Pat_Pathway=TB_Pat_Pathway_6.0
Sample_Size_2.0_TB=length(TB_Pat_Pathway$D1_A2)

##################################
#Initial Classification
###################################

Category_prime=array(100,dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))

Category_prime[which((Facility=="municipal/govt hospital")|(Facility=="municipal health post"))]=6
Category_prime[which(Facility=="local chemist or pharmacy")]=5
Category_prime[which((Facility=="private doctor clinic"|Facility=="Informal provider"|Facility=="private hospital/nursing home"|
                        Facility=="clinic run by NGO"|Facility=="trust hospital"|Facility=="Diagnostic centre")&
                       (Qualification=="BHMS"|Qualification=="BAMS"|Qualification=="BUMS"|
                          Qualification=="traditional healer/quacks/informal provider"|Qualification=="Traditional healer/quack/informal provider"|
                          Qualification=="Traditional healer/quack/informal provider"|Qualification=="BEHMS/CCH/DMLT"|
                          Qualification=="MS Ayurveda"|Qualification=="Not applicable"))]=1
Category_prime[which((Facility=="private doctor clinic"|Facility=="Informal provider"|Facility=="private hospital/nursing home"|
                        Facility=="clinic run by NGO"|Facility=="trust hospital"|Facility=="Diagnostic centre")&
                       (Qualification=="MBBS"|Qualification=="MD chest and TB"|Qualification=="MD (Paed)"|Qualification=="MBBS DA"|
                          Qualification=="MD (DNB)"|Qualification=="MD (General Surgeon)"|Qualification=="MD (Medicine)"|
                          Qualification=="MBBS MD(Orthopedics)"|Qualification=="MD (DCH)"|Qualification=="MD DCH"|Qualification=="MS (ENT)"|
                          Qualification=="M.D.(Path.)"))]=2

#Unknown qualification
Category_prime[which((Facility=="private doctor clinic"|Facility=="Informal provider"|Facility=="clinic run by NGO"|
                      Facility=="private hospital/nursing home"|
                      Facility=="trust hospital"|Facility=="Diagnostic centre")&(Qualification=="1000"))]=3

Categorynew=array(Category_prime, dim=c(Sample_Size_1.0_TB,Stages_Estimation_1.0_TB))

#################################################
#Extraction of Information
#################################################

####################################
#Patient Delay
####################################
TB_Patient_Delay=as.numeric(TB_Pat_Pathway$sym_care)

####################################
#Indicator Variables
####################################

#Test was ordered
Indicator_Provider_Ordered_Test      =array(c(TB_Pat_Pathway$D1_B1,TB_Pat_Pathway$D2_B1,TB_Pat_Pathway$D3_B1,
                                              TB_Pat_Pathway$D4_B1,TB_Pat_Pathway$D5_B1,TB_Pat_Pathway$D6_B1,
                                              TB_Pat_Pathway$D7_B1),dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
#Diagnosis given by the provider at a stage
Indicator_Provider_Diagnosis      =array(c(TB_Pat_Pathway$D1_E1,TB_Pat_Pathway$D2_E1,TB_Pat_Pathway$D3_E1,
                                           TB_Pat_Pathway$D4_E1,TB_Pat_Pathway$D5_E1,TB_Pat_Pathway$D6_E1,
                                           TB_Pat_Pathway$D7_E1),dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
#Switch from the Provider
Indicator_Provider_Switch      =array(c(TB_Pat_Pathway$D1_G7,TB_Pat_Pathway$D2_G7,TB_Pat_Pathway$D3_G7,
                                        TB_Pat_Pathway$D4_G7,TB_Pat_Pathway$D5_G7,TB_Pat_Pathway$D6_G7,
                                        TB_Pat_Pathway$D7_G7),dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
#Treatment given by the provider at a stage
Indicator_Provider_Treatment      =array(c(TB_Pat_Pathway$D1_F1,TB_Pat_Pathway$D2_F1,TB_Pat_Pathway$D3_F1,
                                           TB_Pat_Pathway$D4_F1,TB_Pat_Pathway$D5_F1,TB_Pat_Pathway$D6_F1,
                                           TB_Pat_Pathway$D7_F1),dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))

#######################################
#Tests ordered, performed and Diagnosis given out
#######################################

Tests_Ordered      =array(c(c(TB_Pat_Pathway$D1_B2.1,TB_Pat_Pathway$D1_B2.2,TB_Pat_Pathway$D1_B2.3,
                              TB_Pat_Pathway$D1_B2.4,TB_Pat_Pathway$D1_B2.5,TB_Pat_Pathway$D1_B2.6),
                            c(TB_Pat_Pathway$D2_B2.1,TB_Pat_Pathway$D2_B2.2,TB_Pat_Pathway$D2_B2.3,
                              TB_Pat_Pathway$D2_B2.4,TB_Pat_Pathway$D2_B2.5,TB_Pat_Pathway$D2_B2.6),
                            c(TB_Pat_Pathway$D3_B2.1,TB_Pat_Pathway$D3_B2.2,TB_Pat_Pathway$D3_B2.3,
                              TB_Pat_Pathway$D3_B2.4,TB_Pat_Pathway$D3_B2.5,TB_Pat_Pathway$D3_B2.6),
                            c(TB_Pat_Pathway$D4_B2.1,TB_Pat_Pathway$D4_B2.2,TB_Pat_Pathway$D4_B2.3,
                              TB_Pat_Pathway$D4_B2.4,TB_Pat_Pathway$D4_B2.5,TB_Pat_Pathway$D4_B2.6)),
                             dim=c(Sample_Size_2.0_TB, 6,Stages_Estimation_1.0_TB))

Tests_Performed=array(c(TB_Pat_Pathway$D1_B4,TB_Pat_Pathway$D2_B4,TB_Pat_Pathway$D3_B4,
                        TB_Pat_Pathway$D4_B4,TB_Pat_Pathway$D5_B4),dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))

Diagnosis_Given_By_Provider=array(c(TB_Pat_Pathway$D1_E2.1,TB_Pat_Pathway$D2_E2.1,TB_Pat_Pathway$D3_E2.1,
                                    TB_Pat_Pathway$D4_E2.1),
                                  dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Diagnosis_Given_By_Provider_2=array(c(TB_Pat_Pathway$D1_E2.2,TB_Pat_Pathway$D2_E2.2,
                                      TB_Pat_Pathway$D3_E2.2,TB_Pat_Pathway$D4_E2.2),
                                    dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))

#######################################
# Time durations
#######################################

Diagnostic_duration=array(as.numeric(c(TB_Pat_Pathway$D1_E4,TB_Pat_Pathway$D2_E4,TB_Pat_Pathway$D3_E4,
                                       TB_Pat_Pathway$D4_E4)),dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Duration_with_Provider=array(as.numeric(c(TB_Pat_Pathway$dur_prov1,TB_Pat_Pathway$dur_prov2,
                                          TB_Pat_Pathway$dur_prov3,TB_Pat_Pathway$dur_prov4)),
                             dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Waiting_time=array(as.numeric(c(TB_Pat_Pathway$prov1_prov2,TB_Pat_Pathway$prov2_prov3,
                                TB_Pat_Pathway$prov3_prov4)),
                                      dim=c(Sample_Size_2.0_TB,(Stages_Estimation_1.0_TB-1)))
Treatment_duration=array(as.numeric(c(TB_Pat_Pathway$D1_E5,TB_Pat_Pathway$D2_E5,
                                      TB_Pat_Pathway$D3_E5,TB_Pat_Pathway$D4_E5)),
                         dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))


#############################################
#Estimates on duration with the providers
############################################
Time_Category=array(dim=c(Sample_Size_2.0_TB, Stages_Estimation_1.0_TB))
Delta_sam_correct=array(dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Delta_diffe_correct=array(dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Delta_confirm_correct=array(dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Delta_sam_wrong=array(dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Delta_diffe_wrong=array(dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))
Delta_confirm_wrong=array(dim=c(Sample_Size_2.0_TB,Stages_Estimation_1.0_TB))

#Different scenarios in the TB care pathway

for(i in 1:(Sample_Size_2.0_TB))
{
  k=1
  repeat
  {
    if(Indicator_Provider_Diagnosis[i,k]=="Yes")
    {
      if(Diagnosis_Given_By_Provider[i,k]=="TB"|Diagnosis_Given_By_Provider[i,k]=="MDR TB")
      {
        Time_Category[i,k]=min(Observed(Duration_with_Provider[i,k]),Observed(Diagnostic_duration[i,k]))
        
        if(Indicator_Provider_Switch[i,k]=="Yes")
        {
          Delta_diffe_correct[i,k]=(Duration_with_Provider[i,k])-(Diagnostic_duration[i,k])
          Delta_confirm_correct[i,k]=((as.numeric(TB_Pat_Pathway$diag_tmt))[i])-(Delta_diffe_correct[i,k])
        }
        else
        {
          Delta_sam_correct[i,k]=Treatment_duration[i,k]
        }
      }
      else
      {
        Time_Category[i,k]=min(Observed(Duration_with_Provider[i,k]),Observed(Diagnostic_duration[i,k]))
        
        Delta_diffe_wrong[i,k]=(Duration_with_Provider[i,k])-(Diagnostic_duration[i,k])
        Delta_confirm_wrong[i,k]=((as.numeric(TB_Pat_Pathway$diag_tmt))[i])-(Delta_diffe_wrong[i,k])
      }
    }
    else
    {
      Time_Category[i,k]=Duration_with_Provider[i,k]
    }
    if(k==(Stages_Estimation_1.0_TB)|Diagnosis_Given_By_Provider[i,k]=="TB"|Diagnosis_Given_By_Provider[i,k]=="MDR TB"){break}
    k=k+1
  }
}

#extraction of the existent data points

Delta_diff_correct=Observed(Delta_diffe_correct)
Delta_confirmation_correct=Observed(Delta_confirm_correct)
Delta_same_correct=Observed(Delta_sam_correct)
Delta_diff_wrong=Observed(Delta_diffe_wrong)
Delta_confirmation_wrong=Observed(Delta_confirm_wrong)
Delta_same_wrong=Observed(Delta_sam_wrong)
Twaiting=Observed(Waiting_time)

proc.time() - ptm