COMPARATIVE_ANALYSIS.PY
Contact: sielynanderson@gmail.com
Last Updated 1/20/2015

Table of Contents
1. INTRODUCTION
2. FILES
3. DEPENDENCIES
4. PARAMETERS
5. FUNCTIONS AND OUTPUT
 
1. INTRODUCTION

The comparative analysis script is a Python 3.x script designed to run basic statistical tests and generate plots when given files containing metagenomic abundance data and metadata. It is run from the command line with a single argument which is the filename containing the parameters. Example:

python comparative_analysis.py parameters.sh

It is currently capable of running tests on data in two classes (i.e. control versus experimental). It assumes in the metadata file the experimental samples are denoted by a 1 and the control is marked by a 0. It also assumes all files are tab separated. 

2. FILES

comparative_analysis.py - the main script.

mgprofile.py - class that holds the metadata and abundance data. 

ca_functions.py - contains the functions that perform the computation for the main script. 

Optional Files: 

test_comparative_analysis.py - used to test the comparative analysis script. 
 
testparams.sh - gives an example of how the parameters should be formatted

KO2Sample_T2D.tab and Sample2Class_T2D.tab - real abundance data and metadata respectively, that can be used to show the output of the script, specified in testparams.sh. From Qin et al. 2012 (http://www.ncbi.nlm.nih.gov/pubmed/23023125)

3. DEPENDENCIES

pandas
http://pandas.pydata.org/

scipy.stats
http://docs.scipy.org/doc/scipy/reference/stats.html

sklearn.metrics.pairwise.pairwise_distances
http://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.pairwise_distances.html

matplotlib.pyplot
http://matplotlib.org/api/pyplot_api.html

numpy
http://www.numpy.org/

4. FUNCTIONS AND OUTPUT

Normalization: The only normalization implemented currently is relative. For each sample, the abundances of each attribute are summed. Each attribute abundance is then divided by the total sum so that when all attributes are summed now they add up to 1. 

Enrichment: Performs either a student’s t test or the Wilcoxon rank-sum test on the abundance data. Output is stored in a file p-values.tab which for each attribute contains the attribute name, the p-value, and the directionality (either experimental or control currently). The file is organized from smallest p-value to largest. If an attribute is completely absent in both classes of samples, “n/a” is listed by the attribute name instead of a value. 

PCA: Performs a principal component analysis on the two classes of samples over their abundance data, with blue circles denoting control samples and red triangles for experimental samples. Output is saved in a file pca.png.

PCoA: Performs a principal coordinate analysis on the two classes of samples over their abundance data. Several types of distances are available, all from sklearnmetrics.pairwise.pairwise_distances, including "braycurtis" (default), "jaccard", "manhattan", and "cosine." Output is saved in a file pcoa_(distance type).png

5. PARAMETERS

What follows is a list of the possible parameters that can be passed to the script from the parameters file. 

Keyword: Name
Description: The name of the lab or individual running the script.
Options: n/a 
Required: 

Keyword: Collaborator 
Description: Name of the collaborator. 
Options: n/a 
Required: No

Keyword: Sequence_Type
Description: The type of sequencing used to generate the data. 
Options: n/a
Required: No

Keyword: Year
Description: The year the data was generated. 
Options: n/a
Required: No

Keyword: Abundance_Data
Description: Path to the file containing the abundance data. Assumes file is tab-separated and that the attributes are rows and the samples are columns.
Options: n/a
Required: Yes

Keyword: Metadata
Description: Path to the file containing the metadata. Assumes no headers, that the file is tab-separated, that the experimental samples are denoted by a 1 and the control samples are denoted by a 0, and sample names are rows. 
Options: n/a
Required: Yes

Keyword: Enrichment 
Description: Indiciates whether the enrichment test should be performed to generate p-values. 
Options: Yes, No (default)
Required: No

Keyword: Enrichment_Test
Description: Indicates the type of enrichment test to be performed. 
Options: "T-test" for student's t-test (default), "ranksums" for Wilcoxon rank-sum test
Required: No 

Keyword: Normalization
Description: Indicates whether data should be normalized.
Options: None (default) or Relative
Required: No

Keyword: Output_Dir
Description: Name of the output directory where the files generated should be saved. Defaults to the current working directory. 
Options: n/a
Required: No

Keyword: PCA 
Description: Indiciates whether a principal component analysis should be performed.
Options: Yes, No (default)
Required: No

Keyword: PCoA
Description: Indicates whether a principal coordinate analysis should be performed. 
Options: Yes, No (default)
Required: No

Keyword: Dist_type
Description: If a PCoA is to be performed, indicates what distance metric to be used. 
Options: euclidean, manhattan, cosine, braycurtis (default), jaccard, correlation 
full list here http://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise.pairwise_distances.html 
Required: No