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
1. Input file: RE-GOA takes a `.bed` file input, with which each line shoud have at least 3 colums: `chrom`，`chromStart`，`chromEnd`. Few line of an example input file `exampleinput.bed` are shown as follows:
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
4. Output files in `resultpath` fold includes four files, `exampleinput_BP.txt`,`exampleinput_CC.txt`,`exampleinput_MF.txt`,`exampleinput_genes.txt`. 

`exampleinput_BP.txt`,`exampleinput_CC.txt`,`exampleinput_MF.txt` list out enriched terms in BP, CC, MF, and have 9 colunms, which are described as follows:
```
1.  GO id
2.  GO term name
3.  Num of peaks annotated with the term
4.  Num of peaks located in defined REs
5.  background probability of the term
6.  raw p-value
7.  -log p-value
8.  B-H corrected Q-value
9.  -log Q-value
```

`exampleinput_genes.txt` lists out enriched genes, and have 2 columns, which are described as follows:
```
1.  gene
2.  Num of peaks located in a RE which regulates the gene
```

## Requirements
Python3

Python package: pickle, math, and scipy

It usually takes about 1 minite for computing a `.bed` file and and write all text output files.

## Codes for generate Regulatory Elements Gene Ontology Annotations
codes for training andannotating are availible at `generate RE-GOA`, and datas are availible at https://drive.google.com/file/d/17nYBAGKc2ZZK06mY-9RQPZVSWgAoHxgu/view?usp=sharing

## Citation

If you use RE-GOA or RE-GOA associated resources, please cite

Yurun, et al. Annotating regulatory elements by heterogeneous network embedding. 2021.
