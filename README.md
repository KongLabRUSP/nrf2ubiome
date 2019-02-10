##  Project: Nrf2 BL6 PEITC 16S microbiome data analysis
### Study ID: 
### Scientist: Ran Yin, Renyi Wu
### Data Analysis: Davit Sargsyan 
### Created: 11/28/2018

---    

## Table of Contents
[Daily Logs](#logs)  
[Results](#results)   
[Files](#files)
[References](#ref)   

## Daily Logs<a name="logs"></a>
### 02/09/2019
* Corrected mean plots; added weighted least squares.

### 02/07/2019
* Added interactive plots (plotly) for relative abundances at different tax ranks

### 12/17/2018
* Finished verison 1; rerunning alignment with trimming

### 12/03/2018
* Changed virtual memory size on CVI computer #3 to 64Gb (see Reference #2 below for instructions)

### 12/01/2018
* Initial FastQ processing with DADA2

### 11/28/2018
* 16S FastQ files downloaded to:    
`/datastorage/FastQ_2018/16s/November`    
      
Original files were downloaded form [here](https://www.dropbox.com/sh/5hpqzgdrgqy4y9f/AADP5z39Hl1oi9-L8JLExP46a?dl=0)

## Files<a name="files"></a>
1. ***fastq*** folder contains 80 FastQ files (40 pair-ended samples).    
2. ***docs/legends_16s_11-28-2018.xlsx***: 16S sample legend.    
3. ***docs/Gmail - Fwd_ 16s rRNA sample ID legend from Kong lab.pdf***: email with the link to teh Dropbox folder with all 80 FastQ files.    
4. ***source/nrf2ubiome_dada2_v1.R***: FastQ processing with DADA2 pipeline.

## References<a name="ref"></a>
1. [DADA2 Pipeline Tutorial on GitHub, Benjamin Callahan](https://benjjneb.github.io/dada2/tutorial.html)
2. [Change Window 7 virtual memory Size](https://support.lenovo.com/us/en/solutions/HT002951)
3. [Rutgers On-Demand High Power computing](https://ondemand.hpc.rutgers.edu), run RStudio server (specify number of cores)