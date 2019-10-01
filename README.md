# motif_analysis
Python script for motif analysis for Cserhati, 2019

This script is run on a whole genome sequence, and calculates the score for all possible motifs of length n, ordered in lexicographic order. Since there are four bases, there are 4<sup>n</sup> possible motifs. The score of any given motif is S = (O + E)/(O - E), where O is the observed occurrence of a given motif, and E is its expected occurrence based on the base pair background distribution, which is also determined by the script. S then takes a value between -1 and 1 for under-represented and over-represented motifs, respectively.

The script is run the following way:

python motif_analysis.py -i \<inputfile\> -o \<outputfile\> -s \<species name\> -n \<motif length\>
  
where the flags mean the following things:

-i : the input whole genome sequence in fasta format.<br/>
-o : the output statistics file with the number of contigs, the bp distribution, and the motif sequence, the observed and expected occurrences as well as the score value.<br/>
-s : the species name<br/>
-n : the motif length (usually between 5 and 12 bp)<br/>

A sample of the output file is as follows:

#Species        No. chr.        Genome length   A%      C%      G%      T%<br/>
#GCF_002217835.1_Dobs_1.0_genomic.fna   1935    181868570       0.278804967467  0.221933285975  0.222036923996  0.277224822562<br/>
#Motif  Observed        Expected        Score<br/>
AAAAAAA 257119  23815.5822242   0.83045460594<br/>
AAAAAAC 57298   18957.5905639   0.502788177924<br/>
AAAAAAG 51056   18966.4433466   0.458275306026<br/>
AAAAAAT 88520   23680.6058956   0.577888092376<br/>
AAAAACA 63707   18957.5905639   0.541337096463<br/>
AAAAACC 29016   15090.5502374   0.315722941096<br/>
AAAAACG 23060   15097.5971964   0.208671493717<br/>
