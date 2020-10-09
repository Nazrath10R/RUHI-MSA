#!/bin/bash
#$ -cwd ##setting the workign directory to the directory of this script
cat ./logo.txt
####  Help page  ####
display_usage() {
  echo
  echo -e "Run: analysis.sh [File] [TOLERANCE] [OUTPUT] [THREAD]\n"
  echo
  echo -e "File = psm output from MSFragger. Leave .tsv out"
  echo -e "TOLERANCE  = define which tolerance you want"
  echo -e "OUTPUT  = name of the output file"
  echo -e "THREAD  = number of THREADs to use"
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
  echo
  echo "Enter number of thread or CPUs you want to run RUHI-MSA:"
  read THREAD
  #echo

#### Normal mode ####
else
  PSM=$1
  TOLERANCE=${VARIABLE:10}
  OUTPUT=${VARIABLE:$PSM}
  THREAD=${VARIABLE:2}
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

#Slpit file in THREAD number of files


number_of_lines=$(< "./input/$PSM.txt" wc -l)
echo $number_of_lines
((lines_per_file = ($number_of_lines + $THREAD - 1) / $THREAD))
echo $lines_per_file
gsplit -l $lines_per_file -d "./input/$PSM.txt" $PSM


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
# 3) run RUHI-MSA on parallel on the diferent files
echo "Analysing " $PSM "..."
for each in ./$PSM* ; 
do ./src/RUHI-MSA.run $each $TOLERANCE &
done
wait # wait for all files to run
rm ./$PSM*

#combine files and keep only one header
awk 'NR==1; FNR==1{next} 1' ./Results/"$PSM"*_out.txt > ./Results/"$PSM"_OUT.txt

python3 ./scripts/concatenatefiles.py $PSM

#delete temporary files
for each in ./Results/"$PSM"*"_out.txt" ;
do  rm $each
done

rm ./Results/"$PSM"_OUT.txt


echo "RUHI-MSA is concluded"
echo
echo $PSM"_OUTPUT.txt has been saved in Results"
echo
echo "-------------------------------------------"
echo
echo "Starting plot creation"
echo "..."
# 3) create plots
##PTMs results


##Mass shifts
{
    python3 ./scripts/mass_shifts_plot.py $PSM
} || {
   # command which should be run instead of the above failing command
   python ./scripts/mass_shifts_plot.py $PSM
}
echo
echo "Plot with mass shifts present in the sample saved in:"
echo " ./Results/Myplots/"$PSM"_Mass_Shifts"
echo
echo "..."
echo
{
    # command which may fail and give an error
    python3 ./scripts/PTMs_plot.py $PSM
} || {
   # command which should be run instead of the above failing command
   python ./scripts/PTMs_plot.py $PSM
}

echo
echo "Plot with PTMs present in:"
echo " ./Results/Myplots/"$PSM"_PTMs_PROFILE.html"
