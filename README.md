# Sewer

Sewer pipeline Tutorial

First ensure that you are running as sudo and the required conda enviorment -

1) sudo -s

2) conda activate sewer


after that, you should import the relevent files from data2 to the sewer folder in data3 by the importEnvs script.

make sure that the name in the script args is identical to the folder name in data2

3) ./importEnvs NGSXX_XXXX2021   i.e - NGS84_21042021.

cd to the that new folder that you imported -

4) cd /data3/sewer/NGSXX_XXXX2021

open a new screen - 

5) screen -S NGSXXsewer

run the script -

6) python /data/projects/Dana/scripts/Sewer/pipeline.py BAM/ 5 refs/REF_NC_045512.2.fasta

The parameters are - 
  1] the path for the script
  
  2] the path for the bam folder
  
  3] the minimum coverage to get value
  
  4] path for the reference sequence
  
The results will create in results folder when the script is finished.

The monitored_mutations.csv need to be joined to total_monitored_mutations.csv file, that's the base of the compressed_table.csv

the survellience_table.csv is for the EnvSurf file (Itay & Merav)
