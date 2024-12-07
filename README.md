# Sketching_project


Storing and Querying for Large Collections of Reads:
The code for this component of the project is located in largecollections.py. Make sure that the 5 fasta files are downloaded and stored in the same folder as the code.
Prior to running the code, set up a conda environment with the following installations.<br>
conda create -n bioenv python=3.12<br>
conda activate bioenv<br>
conda install -c conda-forge sourmash-minimal<br>
conda install mmh3<br>
conda install -c bioconda -c conda-forge sourmash<br>
python3 largecollections.py <br>

Benchmarking Bloom, Cuckoo, and Quotient Filters:
The code for this component of the project is located in occupancy_benchmarking.ipynb.
To run this file, you can use the previously configured bioenv, then install pyprobables and matplotlib:
conda install pyprobables<br>
conda install matplotlib<br>

In the main() function, there is a function for each approach called 1) brute_force_storage on line 305, 2) aggregate_bloom_filter on line 309 and 3) sequence_bloom_tree on line 313. Run each approach individually by commenting out the function calls for the other appraoches.

The command line output from the print statements that result from running each appraoch are provided in the following txt files: brute_force_method_output.txt, aggregate_bloom_filter_output.txt, and sbt_output.txt

Run both benchmarking.py (for bloom and cuckoo) and quotient_benchmarking.py (for quotient) to see how we got the simulated data for our implementations. quotient.py was the first implementation of our custom quotient filter. The amount of inserted elements and capacity can be modified to see the trends that appear as occupancy increases. Make sure you install the proper libraries.
conda install mmh3 
conda install memory-profiler

Contributions.txt includes a breakdown of who completed what in our team.

Code Credits: <br>
Adapted from Sourmash Python API for SBT function - https://sourmash.readthedocs.io/en/v2.3.1/api.html<br>
Adapted from UMD JellyFish documentation for K-mer counting - https://www.genome.umd.edu/jellyfish.html<br>
Utilized ChatGPT to assist in creating the SearchFunction class for the SBT function - 
OpenAI. (2023). ChatGPT (Mar 13 version) [Large language model]. https://chat.openai.com/chat


