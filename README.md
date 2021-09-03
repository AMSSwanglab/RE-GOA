# RE-GOA
A tool for non-coding genomic regions function enrichment analysis based on  Regulatory Elements Gene Ontology  Annotation (RE-GOA)
## Installation
1. Run the following commands for installation:
```
wget https://github.com/AMSSwanglab/RE-GOA/archive/master.zip  
unzip master.zip
cd RE-GOA-master
```  
2. Download the necessary files for SpecVar into RE-GOA-master at: https://drive.google.com/file/d/17nYBAGKc2ZZK06mY-9RQPZVSWgAoHxgu/view?usp=sharing and run the fowllowing command:
```
tar -zxvf REGOAd_Data.tar.gz
```  
## Run RE-GOA
1. Inputfile: exampleinput.bed
```
chr1	3389104	3389281
chr1	3991781	3992030
chr1	4333622	4333830
chr1	4402701	4402873
chr1	4424781	4425411
chr1	4463463	4465117
chr1	4469419	4470362
chr1	4527428	4529092
```
2. Edit`configures.txt`, including `path` for path of input file; `datapath=./datamgi/` for mm9 bed file and `datapath=./datahg/` for hg19 bed file; `resultpath` for path of output files; `filename` for input file name
```
path	./inputfiledict/
datapath	./datamgi/
resultpath	./outputfiledict/
filename	exampleinput.bed
```
3. After preparing the input files as above, run the following command:
```
python peaksanalysis.py
```
4. Output files in `resultpath` fold:
