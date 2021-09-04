# RE-GOA
A tool for non-coding genomic regions functional enrichment analysis based on  Regulatory Elements Gene Ontology  Annotation (RE-GOA). Resources for RE-GOA in mouse and human are availible.
## Installation
1. Run the following commands for installation:
```
wget https://github.com/AMSSwanglab/RE-GOA/archive/master.zip  
unzip master.zip
cd RE-GOA-master
```  
2. Download the necessary files for RE-GOA into RE-GOA-master at: https://drive.google.com/file/d/17nYBAGKc2ZZK06mY-9RQPZVSWgAoHxgu/view?usp=sharing and run the fowllowing command:
```
tar -zxvf REGOA_Data.tar.gz
```  
## Run RE-GOA
1. Prepare input file: RE-GOA takes a `.bed` file as input, with which each line shoud have at least 3 columns: `chrom`，`chromStart`，`chromEnd`. `.bed` file for `mm9` or `hg19` are accepted. Few lines of an example input file `exampleinput.bed` are shown as follows:
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
2. Edit`configures.txt`, including:
*  `path` for path of input file; 
*  `datapath=./datamgi/` for `mm9` bed file and `datapath=./datahg/` for `hg19` bed file;
*  `resultpath` for path of output files; 
*  `filename` for input file name.

For example:
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

`exampleinput_BP.txt`,`exampleinput_CC.txt`,`exampleinput_MF.txt` list out enriched terms in Biological Process, Cellular Component, Molecular Function, and have 9 colunms in each line, which are described as follows:
```
1.  GO id
2.  GO term name
3.  m (Num of peaks annotated with the term)
4.  N (Num of peaks located in defined REs)
5.  p (background probability of the term)
6.  raw p-value=Binom(m,N,p)
7.  -log p-value
8.  B-H corrected Q-value
9.  -log Q-value
```
Terms are sorted by `p-value`.

`exampleinput_genes.txt` lists out enriched genes, and have 2 columns, which are described as follows:
```
1.  gene
2.  n (Num of peaks located in a RE which regulates the gene)
```
Genes are sorted by `n`.

## Requirements
Python3

Python package: pickle, math, and scipy

It usually takes about 1 minite for computing a `.bed` file and and write all text output files.

## Codes for generate RE-GOA
Codes for training andannotating are availible at `generate RE-GOA`, and associated datas are availible at https://drive.google.com/file/d/17nYBAGKc2ZZK06mY-9RQPZVSWgAoHxgu/view?usp=sharing

## Citation

If you use RE-GOA or RE-GOA associated resources, please cite

Yurun, et al. Annotating regulatory elements by heterogeneous network embedding. 2021.
