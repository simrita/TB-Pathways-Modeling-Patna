###################################################################################
#This fragment of the code reads the data files, creates the final data-set and then 
#extracts information on the qualification and the facilities of the providers
###################################################################################

# TO clear the memory
rm(list=ls())

ptm <- proc.time()

####################
#Define the parameters
####################

Stages_Data_Entry_TB=7

##############################
#User-Defined Functions
##############################

########################
#To Remove NA entries from an array
########################

Observed=function(x)
{
  x=x[is.na(x)==FALSE]
  return(x)
}

################################
#Reading the Data files
################################

# Original Data Set with 83 patients
TB_Mum_Pathway_1.0 = read.csv("tb_cases_mumbai_28_jan_2015_2015_01_28_11_11_04.csv",colClasses="character")
# Original Data-set of 23 MDR patients (first diagnsis as well as subsequent)  
TB_Mum_Pathway_2.0 = read.csv("mdr_mumbai_28_jan_2015_2015_01_28_11_10_02.csv",colClasses="character")
# Modified Data set of 83 patients with qualification entries revisited
TB_Mum_Pathway_3.0 = read.csv("Pathway to care-FMR provider database with PPIA geo codes.csv",colClasses="character")
# Modified Data set of 76 (TB+ TB-MDR) patients with qualification entries revisited
TB_Mum_Pathway_4.0 = read.csv("Final FMR_provider_mapping_TB_mumbai_provider database_and_patient_geo-coordinates.csv")
# Modified Data set of 23 MDR patients with qualification entries revisited
TB_Mum_Pathway_5.0 = read.csv("Final FMR_provider_mapping_MDR-TB_mumbai_provider database_and_patient_geo-coordinates.csv")

##########################################################################
#Truncating the extra rows
##########################################################################

TB_Mum_Pathway_3.0=subset(TB_Mum_Pathway_3.0, TB_Mum_Pathway_3.0$Patient.ID1!="")
TB_Mum_Pathway_4.0=subset(TB_Mum_Pathway_4.0, TB_Mum_Pathway_4.0$Patient.ID1!="")
TB_Mum_Pathway_5.0=subset(TB_Mum_Pathway_5.0, TB_Mum_Pathway_5.0$Patient.ID1!="")

###############################################
# Merging the data-fields
##############################################

# Data set with 83 patients with all available Qualification entries
TB_Mum_Pathway_6.0=merge(TB_Mum_Pathway_1.0, TB_Mum_Pathway_3.0, by.x="ID2", by.y="Patient.ID2")
#Creating 76 TB and TB-MDR Cases from the above data set with all qualification data points
TB_Mum_Pathway_7.0=merge(TB_Mum_Pathway_4.0, TB_Mum_Pathway_6.0, by.x="Patient.ID2", by.y="ID2")

# #TB Cases
# TB_Mum_Pathway_8.0=subset(TB_Mum_Pathway_7.0, TB_Mum_Pathway_7.0$A6_B=="New pulmonary TB"|TB_Mum_Pathway_7.0$A6_B=="Retreatment")
# TB_Mum_Pathway_9.0=merge(TB_Mum_Pathway_4.0,TB_Mum_Pathway_5.0, by.x="Patient.ID2", by.y="Patient.ID2", all=TRUE)
# #TB+MDR Cases
# TB_Mum_Pathway_10.0=merge(TB_Mum_Pathway_6.0, TB_Mum_Pathway_9.0, by.x="ID2", by.y="Patient.ID2")
# #MDR-cases
# TB_Mum_Pathway_11.0=merge(TB_Mum_Pathway_2.0, TB_Mum_Pathway_3.0, by.x="ID2", by.y="Patient.ID2")

#################################################
#Final Sample Set
#################################################

TB_Mum_Pathway = TB_Mum_Pathway_7.0
Sample_Size_1.0_TB=length(TB_Mum_Pathway$D1_A2)

##################################################
#Extracting the Qualification and Facilities
##################################################

Facility = array(c(TB_Mum_Pathway$D1_A4,TB_Mum_Pathway$D2_A4,TB_Mum_Pathway$D3_A4,TB_Mum_Pathway$D4_A4,
                   TB_Mum_Pathway$D5_A4,TB_Mum_Pathway$D6_A4,TB_Mum_Pathway$D7_A4),dim=c(Sample_Size_1.0_TB,Stages_Data_Entry_TB))

# The two available qualification entries

Qualification_1.0 =  array(c(as.character(TB_Mum_Pathway$D1_A2),TB_Mum_Pathway$D2_A2, TB_Mum_Pathway$D3_A2, 
                                          TB_Mum_Pathway$D4_A2, TB_Mum_Pathway$D5_A2, TB_Mum_Pathway$D6_A2, 
                                          TB_Mum_Pathway$D7_A2),dim=c(Sample_Size_1.0_TB,Stages_Data_Entry_TB))

Qualification_2.0 = array(c(TB_Mum_Pathway$Type.of.facility.1, TB_Mum_Pathway$Type.of.facility.2, 
                            TB_Mum_Pathway$Type.of.facility.3, TB_Mum_Pathway$Type.of.facility.4, 
                            TB_Mum_Pathway$Type.of.facility.5, TB_Mum_Pathway$Type.of.facility.6, 
                            TB_Mum_Pathway$Type.of.facility.7),dim=c(Sample_Size_1.0_TB,Stages_Data_Entry_TB))

#######################################################
# Merging the Qualification Data
#######################################################

Qualification=array(dim=c(Sample_Size_1.0_TB,Stages_Data_Entry_TB))

for(i in 1:Sample_Size_1.0_TB)
{
  k=1
  repeat
  {
    if(Qualification_1.0[i,k]!="Dont know")
    {
      Qualification[i,k]=Qualification_1.0[i,k]
    }
    else
    {
      if(Qualification_2.0[i,k]!="")
      {
        Qualification[i,k]=Qualification_2.0[i,k]
      }
      else
      {
        Qualification[i,k]=1000
      }
    }
    if((k==Stages_Data_Entry_TB)||(Qualification_1.0[i,k+1]=="#NULL!")){break}
    k=k+1
  }
}


proc.time() - ptm

