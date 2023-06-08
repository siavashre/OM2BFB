# OM2BFB
OM2BFB a tool for detecting breakage fusion bridge cycles with OGM technology
## About
Breakage Fusion Bridge (BFB) is a chromosomal aberration occurring during cell division where a chromosome breaks and fuses with another resulting in a bridge-like structure. As the cell divides, the bridge stretches and breaks, leading to abnormal configurations when the fragments rejoin. BFB plays a significant role in genomic instability, and its study is crucial for understanding evolution and disease.
OM2BFB is a tool for detecting BFB events using Optical genome mapping technology. 


## 1. Prerequisites: 
- python3.
- BFBFinder tool. For instructions see [BFBFinder installation](https://github.com/shay-zakov/BFBFinder). 
-  R (version >="4.3")
- `DNAcopy` (R library):  For instructions see [DNAcopy installation](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)
- `numpy` (python library): `pip install numpy`
- `matplotlib` (python library): `pip install matplotlib`
- `pandas` (python library): `pip install pandas`
## 2. Installation:
`git clone https://github.com/siavashre/OM2BFB.git`
## 3. Usage:
`runOM2BFB.py` takes a OM alignments and SVs and CNVs as an input. Here is a explanation:
-  `-r` Path to RefAligner CNV call file (rmcap format)
-  `-c` Path to file contating centromere regions
-  `-n` Project name
-  `-o` Path to folder for saving outputs
-  `-s` Path to RefAligner SV call file (smap format).
-  `-f` Path to FaNDOM SV call. This flag is only reqiured if smap file is not available. 
-  `-x` Path to RefAligner alignment file. (xmap format).
-  `-fol` Path to folder containing molecules to contigs alignment. 
-  `-cmap` Path to bionano contigs file. (cmap format).
-  `-cov` Bionano Sample coverage, default is 77.
-  `-bfbfinder` Path to BFBFinder jar file. 
As an example:
```
python runOM2BFB.py -r snu16.rcmap -c hg38.centro -n SNU16 -o output/snu16 (-s SNU16.smap -f snu16_sv.txt) -x snu16.xmap -fol snu16_alignment/ -cmap snu16.cmap -cov 100 -bfbfinder /Path/To/BFBFinder.jar
```
