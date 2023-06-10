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
-  `-r` Path to RefAligner CNV call file (rmcap format). this file contains Copy number information.
-  `-c` Path to file contating centromere regions. if your reference genome is hg38 please use hg38_centro.txt and if it is hg19 please use hg19_centro.txt
-  `-n` Project name
-  `-o` Path to folder for saving outputs
-  `-s` Path to RefAligner SV call file (smap format). this file contains SV information.
-  `-f` Path to FaNDOM SV call. This flag is only reqiured if smap file is not available. 
-  `-x` Path to RefAligner alignment file. (xmap format). this file contains alignment information.
-  `-fol` Path to folder containing molecules to contigs alignment. OM2BFB uses this for calculating foldback reads multiplicity. In BionanoSolve pipeline it can be find in /contigs/annotation/refine1_ExperimentLabel/
-  `-cmap` Path to bionano contigs file. (cmap format). This is assembled con
-  `-cov` Bionano Sample coverage, default is 77.
-  `-bfbfinder` Path to BFBFinder jar file. It is in the "build/libs/BFBFinder.jar" dir. 
As an example you can download test files for HCC827 and if you install prereqiostes correctly it should take ~5 min to run.
```
python runOM2BFB.py -r test_files/cnv_rcmap_exp.txt -c hg38.centro -n HCC827 -o test_files/output -s test_files/exp_refineFinal1_merged_filter_inversions.smap -x test_files/exp_refineFinal1_merged.xmap -fol test_files/refine1_ExperimentLabel/ -cmap test_files/exp_refineFinal1_merged_q.cmap -cov 100 -bfbfinder /Path/To/BFBFinder.jar
```
## 4. Outputs:
#### ****`amplicon[number]_chr[number]_ans.txt`**** 
This file contains segments information.
- `{Segment}`: Segment name
- `{Chromosome}`: Chromosome number
- {StartPos}: Segment starting position
- {EndPos}: Segment ending position
- {CN}: Segment predicted copy number
- {RightFoldIds}: List of contig IDs that supports right foldback for this segment
- {RightCount}: Predicted multiplicity for right foldback.
- {LeftFoldIds}: List of contig IDs that supports left foldback for this segment
- {lefCount}: Predicted multiplicity for left foldback
- 
#### ****`amplicon[number]_chr[number].csv`**** 
This file is a intermediat file for OM2BFB contains CBS segmentation input.

