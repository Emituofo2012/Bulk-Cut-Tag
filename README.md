# Bulk-Cut-Tag
Bulk cut tag analysis pipeline

Usage:
python Cut-tag-Bulk.V0.2.py --samplecfg sample.cfg.mouse --dataType cleanData --conditionCompare condition --species Mmus --IgGsample IgG.txt

Five Parameters are required, include:

First:
--samplecfg,a text file. it is a sample config file with format as:
samplename  dataPath  condition
1A1 dataPath 1A
1A2 dataPath 1A
three columns,First columns is a sample name,Second columns is a data path which contain samples sequence path,Third is a condition(control or treatment).Seperate by tab.

Second:
--dataType,string. Specify if it is cleanData or rawData. If cleanData is specified, the steps of removing adapters and low-quality reads will not be included.

Third：
--conditionCompare,a text file.Specify the conditions for analyzing differential peaks.Seperate by tab.
format(it means:2A as control and 2C as treatment):
Control Treat
2A  2C

Fourth:
--species,string. Specify the species, this parameter involves genome(.fa), annotation files(.gtf) and other related information.The currently supported species include human and mouse.

Fifth：
--IgGsample,a test file.Specify which conditions are IgG controls for which conditions.
Sample  Blank
2A  1A
1A is a IgG control sample for 2A

After running, a series of shell scripts will be generated. The scripts start with ordinal numbers such as "First", "Second", etc. This is also the order in which the shell scripts should be run.

You can use the “Multiple.processing.py” script to submit these tasks in bulk. This script can help automate the submission and execution of multiple scripts to speed up task processing.
