# The Multispecies Coalescent Model: Theory && Praxis   
## I. Modeling Gene Tree Likelihoods

---
_**Gene Tree**_: a phylogeny inferred from a single aligned set of homologous sequences. 

_**Species Tree**_: a phylogeny inferred from either multiple gene trees _or_ multiple sets of aligned homologous sequence data for several species. 

---
>Modeling the probabilities of multiple gene trees allows us to find the species tree that maximizes the probability of observing that set of gene trees. Before MSC approaches were popular, super trees were a common means of generating species trees, which invovles concating several loci together and treating the resulting data-set as a single, partitioned locus. This approach suffers from several drawbacks, most notably **incomplete lineage sorting**. Incomplete lineage sorting occurs when gene trees do not fully correspond with species tree for a set of taxa. Incomplete lineage sorting can occur when speciation events are recent enough that insufficient time has elapsed for a suitable number of differences to accumulate between two taxa. A growing body of evidence suggests that, in some cases, the genetic differentiation between two species may be confined to certain genomic compartments, probably corresponding with differential adaptation to different environmental and ecological conditions. Gene flow can also exacerbate gene-tree species-tree incongruence. The MSC is far more robust to incomplete lineage sorting due to the fact that it intrinsically accounts for potential differences between gene topologies and branch lengths. 

In traditional molecular phylogenetics, maximum likelihood (ML) or Bayesian approaches are used to model the expected number of nucleotide substitutions (the branch lengths) for a proposed tree topology. The proposed tree topology is tweaked with tree editing algorithms (i.e., subtree pruning && regrafting, nearest neighbor interchange, and tree bisection && regrafting) and branch lengths recalculated until the supposed maximum of some distribution is reached. The time-stationary, time-continuous Markov chain models used to determine the likelihood of a gene tree given an aligned set of sequences is only suitable for the production of gene trees. Rather, a different mathematical model is used, termed **the multispecies coalescent model**. Coalescent theory has its roots in population genetics, where it is often used to trace the likehood of two alleles sharing a common ancestor in a population, among other things. Commonly, coalescent approaches employ a single variable, $\theta$ describing the "effective population size" of a given population. In reality, this value represents the product of some mutation rate $\mu$ and the effective population size:

$$
\theta = 4N_e\mu 
$$

Extensions of this model allow $\mu$ to vary in certain contexts. However, for the time being assume that $\mu$ is the same for all species in a given dataset, which is a reasonable assumption for phylogenies of closely related species. 

In the context of phylogenetics, the multispecies coalescent approach (hereafter: MSC) allows to model the likelihood of observing a given gene tree in terms of speciation times (a vector specifying the 'time' until the next speciation event) and population sizes by hindcasting. Each species or independent lineage in a given gene tree is treated as a population. Two populations coalesce if they share a single ancestral population. Coalescence _can_ occur in any ancestral population that can be the ancestor to two given populations. So, in the gene tree (((A,B),C),D) species A and B can coalesce in either population AB,  population ABC, or population ABCD. By contrast, population C can only coalesce in populations ABCD or ABC. 

In order to use MSC to model the likehood of a gene tree, we need a few equations:
* a function to model the likeihood of a gene tree given some estimate of $\theta$ and species ages (these are branch lengths given in coalescent units and represented by the $\tau$), 
* a function to model the likelihood of coalescence occuring in the time interval that the coalescent ancestral population exists, 
* and a function describing the likeihood of observing a gene tree that matches with a given species tree. 

An important aspect of the MSC process is the coalescent rate, which is $2/\theta$. The units of 'time' are also important, here. One time unit is the time required to accumulate a single mutation per site. 

Each lineage or population in a dataset is traced _backwards in time_ until its termination (i.e., coalescence) at sometime $\tau$. The number of lineages that enter, $m$ and leave, $n$, that lineage at it's termination (i.e., the number of lineages entering and leaving the coalescent population going _backwards in time_) are recorded. Accordingly, if $n \geq 1$, the coalescent processes consists of $n$ disjointed subtrees/lineages. 

Let's consider a single coalescent event that takes the total number of lineages in a given tree from $j$ to $j-1$ lineages. This is ultimately a Poisson Point process, for you statistics nerds out there, meaning the waiting time until the next coalescent event, $t_j$ can be modeled as random variable with an exponential variable, whose probability is defined as 

$$
f(t_j)=\frac{j(j-1)}{2}\frac{2}{\theta}\text{exp}\{-\frac{j(j-1)}{2}\frac{2}{\theta}t_j\}
$$
where $j$ is defined as 
$$
j=m, m-1,...,n+1. 
$$. 

Note that this coalescent waiting time is available directly from the gene tree, not a parameter to be modeled. 

A given lineage ends at some time $\tau$, as discussed above. However, coalescence can occur at any point over this time. The probability that coalescence does not occur **by** $\tau$ is the probability that coalescence does not occur **within** the difference between $\tau$ and and the sum of times associated with this population up to the time of coalescence, $\tau$. This time interval is, in short, the time at which the population exists, and is given by
$$
\tau -[\ t_m + t_m-1 + ... + t_{n+1}]\
$$
or, more succinctly, 
$$
\tau - [\ \sum_{j = m}^{n+1} t_j ]\ . 
$$
Given that we're tracing coalescence back in ward in time, we actually _want_ the probability that no coalescence occurs over this time period. This is because the coalescence should occur only at the end of this period, not before it. This probability is given by
$$
\text{exp} \{ -\frac{n(n-1)}{\theta}(\tau - \sum_{j = m}^{n+1} t_j\} . 
$$ 
Thus far, we have functions for deriving estimates for both $\theta$ and $\tau$ values, which we can represent as a single set $\Theta=\{\tau_i , \theta_i\}$. In other words, we've derived $f(t_i|\Theta)$. This is well and good, but we still haven't accounted for the contribution of gene tree topology, which may or may not actually correspond with that of the species tree, which is our ultimate goal here. Once again, let's consider a coalescent event that occurs in $j$ species or lineages. The probability that two lineages in $j$ coalesce is $1/\binom{j}{2}=\frac{2}{j(j-1)}$. Briefly, $\binom{j}{2}$ is the number of potential coalescent events that can occur for $j$ lineages. If each pair of lineages coalesce with same probability, then the probability of any two lineages coalescing is $1/\text{no. of potetial coalescent events}$, or $1/\binom{j}{2}$. 

Multiplying all of these probabilities together gives the probability of observing a given gene tree
$$
\prod_{j=n+1}^{m} [\ -\frac{2}{\theta}\text{exp} \{ \frac{j(j-1)}{\theta}t_j \} ]\ \text{exp} \{ -\frac{n(n-1)}{\theta} (\tau - \sum_{j=m}^{n+1}t_j)\}.
$$

Taking the product of this across the entire tree gives a function for it's likelihood given some our estimates of $\Theta$.

Given that a set of gene trees is required for inference using MSC, we often describe $G$ as being as set of gene trees: $G=\{G_i\}$, where each gene tree, $G_i$, contains both a topology, $T$ and branch lengths, $t$. Thus, $G_i=\{T_i,t_i\}$. The product of our likelihood function across all gene trees, $G$ in our data set gives: 
$$
f(G|\Theta)=\prod_{i} f(G_i|\Theta)= \prod_{i} f(T_i,t_i|\Theta).
$$

Because of the complicated nature of this algorithm, most programs split up gene trees into quartetts to make calculations easier. 

## II. Maximum Likelihood Analysis: two step methods 
The algorithm that searches species-tree space is otherwise comparable to those used for maximum likelihood inference of gene-trees. Briefly, a starting tree is proposed and parameters are estimated by maximizing the likelihood of observing a gene tree given a species tree topology. The tree is then altered, resulting in a new tree. If the original tree has a lower maximum likelihood, the new tree is accepted and the process is repeated on the new tree. Otherwise, the process is repeated on the same tree until a better tree is found. 

ML methods often require two steps to generate a species tree from a set of aligned loci. Inferring species trees directly from a set of alignments would require integration across all possible gene trees and parameters. This is intractable. As such, ML methods for inferring species trees do so from gene trees. For clarity, a simple work flow (incl. all steps from alignment to species tree) might look something like this:  
1. Align each .fasta file corresponding to an individual locus (using MAFFT, MUSCLE, ClustalW, etc). 
2. Infer a maximum likelihood gene tree for each of the aligned loci (using IQTree2, RAxML-NG, etc). 
3. Feed gene trees into program to infer species tree (ASTRAL is currently the industry standard program for ML based species tree inference).   

Note that in the above, only ML methods are used. **Do not combine Bayesian and Maximum Likelihood approaches, ever**... unless a statistician tells you it's okay first. 

## III. Bayesian Inference: the better option 
It's outside of the scope of this document to cover all of the differences between Bayesian and ML inference. See `ML_Bayesian.md` for a refresher on Bayes' theorem. Briefly, the posterior probability $f(S,G,\Theta|D)$ describes the likelihood of observing species tree $S$, gene trees $G_i=\{T_i, t_i\}$ and population est./coalescent times $\Theta$ given a set of aligned sequence data, $D=\{D_i\}$. Bayes' theorem implies that this posterior distribution is proportional to the product of the likelihood function, $f(D|G)f(G|S, \Theta)$, and posterior probabilities $f(\Theta)f(S)$. These posterior probabilities are distributions describing our _a priori_ expectations about the average value and variance of the estimated parameters. $f(S)$ is unique in that it represents a prior on the topology of the species tree. Often times, this is treated as a Yule model, which fits one parameter for speciation rate. A natural departure from the Yule prior is the Birth-Death model, which fits an additional parameter for the rate of species decay, making it more accurate under certain conditions. 

In short, the Markov Chain Monte Carlo (MCMC) algorithm used in Bayesian analysis roughly 'maps' the posterior distribution by calculating it's value for various estimates. If the estimates increase the value, that state is saved and record. Otherwise, the proposed change of state for the chain (which you can think of as the position on the posterior distribution) is rejected. 

Due to the utility of the MCMC algorithm, species trees can be inferred directly from aligned loci. Moreover, the flexible nature of the MCMC algorithm allows for a number of additional analyses that require more (sometimes many more) parameters. **Because Bayesian methods allow for inference of species trees directly from alignments, they are expected to be more accurate**. Valuable information along the transition from alignment to gene tree in maximum likelihood approaches. As such, ML approaches will require MORE loci to infer an accurate species tree. Moreover, because the gene and species are simultaneously coestimated in Bayesian approaches, parameter estimates are those that optimize gene and species trees, resulting in a higher accuracy for both. 