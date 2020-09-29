#!/bin/bash
#$ -cwd ##setting the workign directory to the directory of this script
cat ./logo.txt
####  Help page  ####
display_usage() {
  echo
  echo -e "Run: analysis.sh [File] [TOLERANCE] [OUTPUT]\n"
  echo
  echo -e "File = psm output from MSFragger. Leave .tsv out"
  echo -e "TOLERANCE  = define which tolerance you want"
  echo -e "OUTPUT  = name of the output file"
  echo

}

if [[ ( $1 == "--help") || ($1 == "-h") ]]
then
  display_usage
  exit 0
fi

### if less than three arguments supplied, display usage  ##
if [[ ( $1 != "--help") &&  $1 != "-h" && ( $1 != "--interactive") &&
  $1 != "-i" && $# -le 0 ]]
then
  display_usage
  exit 1
fi


####  Interactive mode  ####

if [[ ( $1 == "--interactive") ||  ($1 == "-i") ]]
then
  echo
  echo "Runing Interactive mode"
  echo
  echo "Enter name of psm.tsv to analyze:"
  read PSM
  echo
  echo "Enter tolerance:"
  read TOLERANCE
  echo
  echo "Enter output name:"
  read OUTPUT
  #echo

#### Normal mode ####
else
  PSM=$1
  TOLERANCE=${VARIABLE:10}
  OUTPUT=${VARIABLE:$PSM}
fi

#setting default value for tolerance
##TOLERANCE=${VARIABLE:=10}

if [[ ( $1 == "--directory") || ($1 == "-d") ]]
then
  echo
fi



# 1) prepare psm.tsv file
echo
echo
echo "Starting preparation of " $PSM "..."
echo "..."
{
    # command which may fail and give an error
    python3 ./scripts/prep_psm.py $PSM
} || {
   # command which should be run instead of the above failing command
   python ./scripts/prep_psm.py
}

echo "Finished preparing " $PSM
echo
echo "-------------------"
echo
# 2) compile RUHI-msa
echo "Initializing compilation..."
g++ -std=c++17 ./src/main_pedro_.cpp
echo "..."
echo "Compilation finished"
echo
echo "-------------------"
echo
# 3) rename file
mv ./a.out ./src/RUHI-MSA.run
# 3) run RUHI-MSA
echo "Analysing " $PSM "..."
./src/RUHI-MSA.run $PSM $TOLERANCE
echo "..."
echo
echo "RUHI-MSA is concluded"
echo
echo $PSM "_output.txt has been saved in Results"
echo
echo "-------------------"
echo
echo "Starting plot creation"
echo "..."
# 3) create plots
##PTMs results


##Mass shifts
{
    # command which may fail and give an error
    python3 ./scripts/mass_shifts_plot.py $PSM
} || {
   # command which should be run instead of the above failing command
   python ./scripts/mass_shifts_plot.py $PSM
}
echo
echo "Plot with mass shifts present in the sample saved in:"
echo " ./Results/Myplots/"$PSM"_Mass_Shifts "
echo
echo "..."
echo
Rscript ./scripts/PTMs_results_plot.R $PSM
echo
echo "Plot with PTMs present in:"
echo " ./Results/Myplots/"$PSM"_PTMs.pdf"
