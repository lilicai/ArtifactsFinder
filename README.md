### Introduction

Before sequencing, artifacts could be maked every step at the experiment, these artifacts were maked by palindromic structure in the genome. The scripts will search artifacts in given regions at different experimental treatments. 

ArtifactsFinderIVS.py can find IVS (inverted repeat sequences) artifacts, ArtifactsFinderPS.py can find PS (palindromic sequences) artifacts.


### Manual

**Dependencies:**
* Required input files：bed file
* Required software：python3 ，pyfasta，pysam，collections，line_profiler
* Required data:Human genome file (hg19 or hg38)
  
  hg19
  https://hgdownload2.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

  hg38
  https://hgdownload2.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

**Install Using Bioconda**

* conda create -n python3 -c conda-forge -y
* conda install python=3
* conda install bioconda::pyfasta
* conda install bioconda::pysam
* conda install conda-forge::line_profiler
* git clone https://github.com/lilicai/ArtifactsFinder


### Useage

**To run ArtifactsFinder**

	python3 ArtifactsFinderIVS.py -g genome -b test.bed

	python3 ArtifactsFinderPS.py -g genome -b test.bed
	

Note : The genome is the hg19 or hg38 sequences. The bed is in 0-base format. If your region is too long, you can use bedReform.py reform it.

**BED File**
```
An example of a BED File can be viewed here: ([`test.bed`](ArtifactsFinder/example/test.bed)).
```

**IVS blacklist**
```
An example of a IVS blacklist can be viewed here: ([`test.bed.out`](ArtifactsFinder/example/test.bed.out)).
```

**PS blacklist**
```
An example of a PS blacklist can be viewed here: ([`out.backlist.txt`](ArtifactsFinder/example/out.backlist.txt)).
```
Note : 

If you want to eliminate the artifacts caused by palindromic structures, you can use ArtifactsFinder to generate a blacklist of noisy sites in regions of interest. This blacklist can be incorporated into the bioinformatics pipeline to filter the results of variant detection, thereby removing the artificial traces in the specified regions.

###Contact
If you have any questions, please contact lilicai@chosenmedtech.com .
