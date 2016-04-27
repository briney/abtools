# AbTools  
  
AbTools is a package of several tools for performing common antibody repertoire analysis functions:  
`AbCorrect` corrects sequencing errors using either molecular barcodes or homology-based clustering  
`AbFinder` identifies sequences from NGS datasets with similarity to known antibody sequences  
`AbPhylogeny` generates publication-quality phylogenetic trees from antibody sequences  
`AbCompare` provides methods for repertoire-level comparison of antibody sequence data  
  
  - Source code: [github.com/briney/abtools](https://github.com/briney/abtools)  
  - Documentation: [abtools.readthedocs.org](http://abtools.readthedocs.org)  
  - Download: [pypi.python.org/pypi/abtools](https://pypi.python.org/pypi/abtools)  
  
### install  
`pip install abtools`  
  
### requirements  
Python 2.7 (3.x probably doesn't work, but hasn't been tested)  
biopython  
celery  
ete2  
matplotlib  
numpy  
pandas  
seaborn  
pymongo  

All of the above dependencies can be installed with pip, and will be installed automatically when installing AbTools with pip. If you're new to Python, a great way to get started is to install the Anaconda Python distribution (https://www.continuum.io/downloads), which includes pip as well as a ton of useful scientific Python packages.  
      
## AbCorrect  
  
To calculate centroid sequences using unique antibody IDs (UAIDs) that have already been parsed by AbStar:  
`abcorrect -d <database_name> -c <collection_name> -t <temp_dir> -o <output_dir>`  
  
Same as above, but using a remote MongoDB server:  
`abcorrect -d <database_name> -c <collection_name> -i <ip> -p <port> -u <username> -p <password> -t <temp_dir> -o <output_dir>`  
  
Calculate consensus sequences using UAIDs (20 nucleotides long) that haven't been pre-parsed with AbStar:  
  
`abcorrect -d <database_name> -c <collection_name> -t <temp_dir> -o <output_dir> --parse-uaids 20 --consensus`  
  
Same as above, but UAID is 8 nucleotides long and is at the end of each sequencing read:
  
`abcorrect -d <database_name> -c <collection_name> -t <temp_dir> -o <output_dir> --parse-uaids -8 --consensus`  
  
Calculate consensus sequences using homology-based clustering, at a threshold of 96% identity:  
  
`abcorrect -d <database_name> -c <collection_name> -t <temp_dir> -o <output_dir> --identity 0.96 --no-uaids --consensus`  
  
  
## AbFinder    
  
To generate an identity/divergence plot using known antibody sequences in a file named 'known.fasta':    
  
`abfinder -d <database_name> -c <collection_name> -t <temp_dir> -o <output_dir> -d <path/to/known.fasta>`  
  
Same as above, but without updating the MongoDB database with identity information:  
  
`abfinder -d <database_name> -c <collection_name> -t <temp_dir> -o <output_dir> -d <path/to/known.fasta> --no-update`  
  
If known.fasta contains nucleotide sequences (default is amino acid):  
  
`abfinder -d <database_name> -c <collection_name> -t <temp_dir> -o <output_dir> -d <path/to/known.fasta> --no-update --nucleotide`  
  
  
## AbCompare    
  
To compare two MongoDB collections of sequences (collection1 and collection2):  
  
`abcompare -d <database_name> -1 collection1 -2 collection2 -o <output_dir>`  
  
To iteratively compare collection1 to every other collection in the database:  
  
`abcompare -d <database_name> -1 collection1 -o <output_dir>`  
  
To iteratively compare collection1 to every other collection in the database that starts with 'abcd':  
  
`abcompare -d <database_name> -1 collection1 --collection-prefix abcd -o <output_dir>`  
  
To compare two collections (collection1 and collection2) using Jensen-Shannon similarity (default is Marisita-Horn):  
  
`abcompare -d <database_name> -1 collection1 -2 collection2 -o <output_dir> --similarity-method jensen-shannon`  
    
  
## AbPhylogeny  
  
Generate a phylogenic tree from a FASTA-formatted file of sequences (input.fasta), a root sequence (root.fasta):  
  
`abphylogeny -i <path/to/input.fasta> -r <path/to/root.fasta> -o <output_dir>`  
  
  
  
