# Sketching_project


Storing and Querying for Large Collections of Reads:
The code for this component of the project is located in largecollections.py. Make sure that the 5 fasta files are downloaded and stored in the same folder as the code.
Prior to running the code, set up a conda environment with the following installations.<br>
conda create -n bioenv python=3.12<br>
conda activate bioenv<br>
conda install -c conda-forge sourmash-minimal<br>
conda install mmh3<br>
conda install -c bioconda -c conda-forge sourmash<br>
python3 largecollections.py =<br>

In the main() function, there is a function for each approach called 1) brute_force_storage on line 305, 2) aggregate_bloom_filter on line 309 and 3) sequence_bloom_tree on line 313. Run each approach individually by commenting out the function calls for the other appraoches.

The command line output from the print statements that result from running each appraoch are provided in the following txt files: brute_force_method_output.txt, aggregate_bloom_filter_output.txt, and sbt_output.txt
