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

Stages_Data_Entry_TB=4

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
TB_Pat_Pathway_1.0 = read.csv("tb_cases_patna_28_jan_2015_2015_01_28_11_00_29.csv",colClasses="character")


#################################################
#Final Sample Set
#################################################

TB_Pat_Pathway = TB_Pat_Pathway_1.0
Sample_Size_1.0_TB=length(TB_Pat_Pathway$D1_A2)

##################################################
#Extracting the Qualification and Facilities
##################################################

Qualification_1.0 = array(c(as.character(TB_Pat_Pathway$D1_A2),TB_Pat_Pathway$D2_A2, TB_Pat_Pathway$D3_A2, TB_Pat_Pathway$D4_A2),dim=c(Sample_Size_1.0_TB,Stages_Data_Entry_TB))
Facility = array(c(TB_Pat_Pathway$D1_A4,TB_Pat_Pathway$D2_A4,TB_Pat_Pathway$D3_A4,TB_Pat_Pathway$D4_A4),dim=c(Sample_Size_1.0_TB,Stages_Data_Entry_TB))

#######################################################
# Merging the Qualification Data
#######################################################

Qualification=array( dim=c(Sample_Size_1.0_TB,Stages_Data_Entry_TB))

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
      
      Qualification[i,k]=1000
      
    }
    if((k==Stages_Data_Entry_TB)||(Qualification_1.0[i,k+1]=="#NULL!")){break}
    k=k+1
  }
}

proc.time() - ptm

