## Bioinformatics Server to analyse protein sequences

#### This server is build using Python programming coupled with its Web Development library - Streamlit and combining biological analysis libraries like Biopython,py3DMol and stMol to focus on several important analytical techniques to analyse proteins and amino acids some of the features which can be computed from this server are:

- Amino Acid FASTA file uploader.
- Calculation of Physio-Chemical properties of protein using EXPASY programming.
- Multiple Sequence Alignment and BLAST on query uploaded sequence to obtain similar sequences for phylogeny
- In built protein structure visualiser through PDB code and separate locator for desired amino acid.
- Protein-Protein or Protein-Chemical interaction using STRING server API to analyse protein connections and Enrichment analysis along with GO terminologies

---

## Dependancies to install : Perform in python/conda virtual environment 

* Biopython 
```sh
pip install biopython
