# RUHI-MSA

Steps to retrieve data from RetroMiner and populate into MongoDB for RTPEA

## RetroMiner

The last script of RetroMiner is "Data_filtering.sh", which uses a custom version of PeptideShakers' ReportCLI to write out tables with the protein identification information we need for all proteins 

![alt text](https://github.com/Nazrath10R/RetroMiner_to_RTPEA/blob/master/images/RetroMiner%20to%20RTPEA.png)


### Softwares 

* ***SearchGUI*** - http://compomics.github.io/projects/searchgui.html
* ***PeptideShaker*** - http://compomics.github.io/projects/peptide-shaker.html

will be downloaded with the installation

## Prerequisites

| **Platform** | **Languages** | **source**  |	**info**	|	
|----------|--------|-------------|--------------|
| Linux    | bash   | [jq](https://stedolan.github.io/jq/) for bash      | 'sudo apt-get install jq' (linux) or 'brew install jq'	|
| Mac      | R      | [ssh_keys](https://www.digitalocean.com/community/tutorials/-to-configure-ssh-key-based-authentication-on-a-linux-server) for HPC |	follow instructions for your own HPC	|



### Get RT data out - create_results_table.sh (working)

<span style="color:blue">Apocrita</span>

```
sh custom_report2.sh
```

this script runs the modified data export for PeptideShaker and outputs .txt files
loop requires a list of PXDs in a text file. write a way to automatically find the ones not parsed. 
(maybe using reanalysis log) 

run second part of the script for the filtration
awk commands filter RT protein lines out super fast and create LINE and HERV .txt files 

### Create output table 

<span style="color:blue">Apocrita</span>

creates one output table with all results

```
Rscript parser_argumented.R --PXD "PXD00xxxx"
```

makes new folder called final
move all filtered files into final folder


### collate output table 

<span style="color:blue">Apocrita</span>


```
Rscript make_output_table.R
```

## visualise results 

<span style="color:blue">Apocrita</span>

```
Rscript result_interpretation.R
```


## Table Data

<span style="color:blue">DropBox</span>

Table data

```
/data/results/
```

### add metadata

<span style="color:blue">DropBox</span>

```
Rscript adding_consequence_to_output_table.R
````

### convert results to json

<span style="color:blue">DropBox</span>

using example.json

```
Rscript convert_results_to_json_working_final.R
```

move all json files to new folder and put this script in there

```
python Fix_Json_ORF.py
```

move them out of the newly generated folder back out and delete generated folder 



## ProtVista data

<span style="color:blue">DropBox</span>

/variants/

using visualise_example_new.json

```
Rscript conversion_4.R
```

need to add a newly compiled PXD list


### merge into one file per protein 

<span style="color:blue">Apocrita</span>

/ProtVista/results/

```
sh combine_protvista_json_files.sh
```

## ideogram data

<span style="color:blue">DropBox</span>

in data/variants/

using sequences from data/variants/sequences

```
Rscript ideogram.R
```


----------------
### 4. Run RetroMiner

<span style="color:blue"><!-- DropBox --></span>

To run RetroMiner the only command needed is:

```
sh retrominer.sh [PXD00xxxx] [ANALYSIS] [THREADS]
```

where ANALYSIS = 1 (frontend5) ; 2 (frontend6) ; 3 (sm11)

and THREADS = number of cpu cores to use

Please be considerate when using shared nodes and use [htop](https://hisham.hm/htop/) to monitor other users and jobs.

If no arguments are input, the usage is displayed in the terminal (or -h is used)
Using -i, loads up the interactive version guiding you through the input arguments

All proteomic outputs are stored in outputs/PXD00xxxx/ and protein tables in reports/PXD00xxxx/

Also find any filtered protein outputs in results/PXD00xxxx/ 

All console output is written into log files (see logs/PXD00xxxx/) 

Please use RetroMiner_to_RTPEA (ref) for downstream conversion into jSON files and database population
as well as creating Rdata tables for easier visualisations.


