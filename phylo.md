# Primer on Phylogenetics: Theory and Praxis 
## Introduction: What is molecular phylogenetics? 
Molecular phylogenetics is the process by which a set of **homologous** DNA or protein sequences are used to infer evolutionary relationships between related taxa. A phylogenetic tree represents these relationships. The "leaves" of the tree represent extant taxa, whereas internal "nodes" of the tree represent the common ancestors of extant taxa. The length of the branches connecting trees to internal nodes represent the genetic difference that accumulates along that branch. 

Phylogenies can be generated through several different methodologies, each of which maximizes a different criterion to find a "correct" tree and uses a different algorithm. Here are some of the most common methods: 
* **Distance based methods**: Distance based phylogenetic inferences uses observed differences in a provided alignment to calculate a simple measure of genetic distance between extant taxa. These are stored in a data matrix which is then used to construct a distance tree. Because these methods only account for observed differences, they do not account for **homoplasy** (i.e., when two taxa both have a particular nucleotide or amino acid at a given position that arose by convergent evolution rather than from a common ancestor) and may underrepresent the divergence between taxa. 
* **Maximum Parsimony methods**: Parsimony based methods find a tree that the minimizes the total number of supposed mutation events. In other words, parsimony softwares operate on Occam's razor: the tree that requires the least ammount of assumptions (here, presumed changes) is the most likely. 
* **Maximum Likelihood methods**: Maximum likelihood based methods find the tree that maximizes the probability of observing the given data. More on these later. 
* **Bayesian methods**: Bayesian finds the tree that has the highest probability in the posterior distribution. The posterior is made up of both the maximum likelihood distribution (see above) and prior information (these are called "priors"). Accordingly, these methods find the maximally likely tree given both sequence data and prior information. 

In practice, inferring a phylogeny requires two steps: 
1. Alignment 
2. Inferring the tree

Below, I provide a tutorial for both of these steps. The outline I've provided above is precursory and certain phylogenetic analyses may require more background information than an exercise of this scope can provide. 

## Alignment: 
1. The alignment process will require two steps: I. aligning the sequences, and II. removing spurious nucleotides, gaps, etc. from the alignment (or Quality Control/QC). Each of these requires a different software package. To install these, you will need Anaconda installed and configured on your computer. Once Anaconda has been installed and configured, you may install the required software. Here, we will be using MAFFT to align sequences and TrimAl to QC the resulting alignment: 
```
conda install -c bioconda mafft 
conda install -c bioconda trimal
```
2. Run the following: `git clone https://github.com/CorbinBryan/PhyloPractice.git`. This command will clone this github repository into a separate subdirectory of your computer. 
3. Next, change into this subdirectory by doing the following: `cd PhyloPractice`
4. In this repository, you should see a fasta file named `SectAmanitaITS.fa`. This file contains the sequences that we will align, and from which we will infer a maximum likelihood phylogeny. The following code uses MAFFT to generate an alignment and stores that alignment in `ITS_al.fa`. 
``` 
mafft --auto SectAmanitaITS.fa > ITS_al.fa 
```
5. Once the alignment process is finished, you will need to trim it with TrimAl. Here is the appropriate code to do so. Note that the trimmed alignment is saved as `ITS_trim_al.fa`.  
```
trimal -in ITS_al.fa -out ITS_trim_al.fa 
```

Now, your alignment is ready for phylogenetic inference. But wait! Don't move on quite yet. It is imperative that you visually examine your alignment to scan for any poorly aligned sequences. These are often pretty easy to spot. Using a visual alignment viewer software (I use AliView, though others are doubtless available). As in all of bioinformatics, garbage in makes garbage out. As such, one should always check the quality of their alignment before proceeding with phylogenetic inference. 

## Phylogenetic Inference: 
For the purposes of this exercise, we will be inferring the maximum likelihood tree for the ITS gene across taxa in sect. _Amanita_ (an infrageneric taxonomic division of genus _Amanita_). Maximum likelihood approaches require a model of sequence evolution to be selected. These models have different numbers of parameters. These parameters are used to determine the expected number of nucleotide differences between two taxa. The number of parameters in each model depends on the number of different mutation processes it accounts for. The simplest model, the Jukes Cantor model (abbrev. JC69) has just one parameter, which represents the mutation rate for any given nucleotide. By contrast, HKY85, another DNA substitution model, differentiates between the expected rates of transition mutations (purine -> purine or pyrimadine -> pyrimadine) and of transversion mutations. The most complicated model, the general time reversible model, fits a rate parameter for each of the four possible nucleotide transitions. There is some nuance to substitution models that I have not covered here, and many modifications have been made such models to account for variation in evolutionary rates among sites. Some substitution models may assume nucleotides have equal frequencies, while others do not.  

The accuracy and utility of a phylogeny is highly contingent on two things: I. appropriate model selection, and II. the quality of the alignment it accepts as input. We've already seen how to generate a high quality alignment. How do we select our model, then? Well, this is one less cut and dry. The most complicated models may 'over-fit' an alignment, meaning they account for more variability than is actually present in the data and may be misleading. Meanwhile, models that are too simple for a given alignment may be 'over dispersed'. In statistics speak, this means that the model accounts for less variability than is actually present in the data, resulting in a phylogeny that does not accurately reflect the true genetic difference between taxa. As such, one should always choose the model that accurately reflects an alignment with as few parameters as possible. Given that there are a large number of substitution models with an even larger number of potential modifications, this is still a Herculian task. Thankfully, model selection tools have been developed for this very purpose. In this exercise, we will use one such tool called ModelFinder. This tool is automatically implemented in the phylogenetics software we will use, IQTree2. 
1. Install IQTree2: 
```
conda install -c bioconda iqtree
```
2. Run IQTree2 with model finder enabled on your alignment using SH-aLRT and ultra-fast bootstrap as a measure of branch support. These measures of branch support effectively provide a measure of confidence in a given clade. You can think of them as a measure of how frequently a certain clade of taxa appears on the tree. Here's the code: 
```
iqtree2 -s ITS_trim_al.fa --alrt 1000 -B 1000
``` 
3. This will provide you with a completed phylogeny. For this exercise, your tree file should be named `ITS_trim_al.fa.treefile`. Open this using FigTree (or R's `ggtree` package, if you'd prefer). 