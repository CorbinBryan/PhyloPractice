## Getting Started
Before proceeding, ensure that you have the following installed: 
* Visual Studio Code
* Docker
For this exercise, we will be running our analysis in a Docker container which I have prepared. Docker containers are isolated systems that allow you to explore a program without giving it access to the rest of your computers files, allowing for programs to be run without dependency conflicts.

Open a terminal session (if on Windows, try Windows power shell) and type the following. 

```
docker run -it multimuscaria/wi_amanita:1.1
```
## Analysis
1. Run `ls` to view the files in the container. You should see a file called `wi_amanita.fasta`. Examine it's contents by typing `head wi_amanita.fasta`. Note that `head` is a command that shows you the first 20 or so lines of a file. You will notice that each sequence is preceeded by a line that begins with the character `>`. The lines that begin with `>` contain the **sequence IDs**. 
2. Note that the sequence IDs all lack spaces. Certain downstream programs may not run if sequence IDs contain spaces. Normally, I would truncate the names of this fasta file. 
> Phylogenetics can only be performed on aligned loci. Nucleotides present in one taxa may not present in the homologous region of a relative due to mutational point processes (i.e., insertions and deletions, colloquially, in-dels). In-dels and nucleotide substutions must be accounted for by models of DNA evolution used to infer phylogenetic relationships. To infer the presence of either, each site in each sequence must correspond to the exact orthologous site in it's relatives. A site not present in a taxon is represented by a gap (typically given by the chacter `-`). Alignment softwares align sequences by inferring the presence and location of gaps. The resulting alignment can be thought of as a matrix: each row represents a taxon and each column represents a different orthologous site in the sequence. Each sequence starts and stops at the same site and is of equal length after alignment. 
3. We'll align our phylogeny with a program called MAFFT. This is a good time to introduce a key concept in bash scripting: flags. Flags are how additional arguments or parameters are passed to bash commands. Short-hand flags are often a single letter. Try the example below 
```
# Example 1: Run this and evaluate the output. 
ls
# Example 2: Notice how the addition of the -l flag modulates output. 
ls -l 
```
4. Flags can also be longer. The `--auto` flag for the program MAFFT is proceeded by two dashes, as is customary for long-form software inputs. Note that MAFFT does not automatically save the outputed alignment file. You will need to do that manually. You can save the outputs of any program that prints to standard out (i.e., prints the output in your terminal window) using the following syntax: `[PROCESS] > [OUTPUT FILE NAME]`. With these concepts in mind, let's align our loci:
```
mafft --auto wi_amanita.fasta > wi_amanita_al.fasta
```
5. In the prior step, we saved our alignment to `wi_amanita_al.fasta`. We will now perform phylogenetic inference on this alignment. The `-s` flag of iqtree2, the program we will be using to generate our phylogeny, tells the program to perform tree search (infering the maximum likelihood phylogeny). By default, the program will also run ModelFinder, which finds the best nucleotide substitution model for a given set of sequences. Run the code below to generate your tree. 
```
iqtree2 -s wi_amanita_al.fasta 
``` 
6. Unfortunately, the tree we just made has no support values, and thus we have no measure of how "probable" our evolutionary groupings are. Rerun your analysis using SH-aLRT as a measure of branch support with 1000 iterations. Do not use boot strapping. To get started run `iqtree2` to bring up the iqtree help screen. Use the information provided to determine how to do this on `wi_amanita_al.fasta`
7. User `docker cp` to copy your the `.treefile` file from the previous step. Determine how to do this using the documentation available online (i.e., google "docker cp documentation"). Let me know if you can't get this part (no shame, docker is hard). 