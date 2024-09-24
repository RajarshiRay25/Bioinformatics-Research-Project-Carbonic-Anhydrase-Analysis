# Bioinformatics-Research-Project-Carbonic-Anhydrase-Analysis
## Bioinformatics Server to analyse protein sequences
---
### Public access link for using our developed server.

* Seq-Eclipse : [Visit our server](https://seqeclipse.streamlit.app/)

#### This server is build using Python programming coupled with its Web Development library - Streamlit and combining biological analysis libraries like Biopython,py3DMol and stMol to focus on several important analytical techniques to analyse proteins and amino acids some of the features which can be computed from this server are:

- Amino Acid FASTA file uploader.
- Calculation of Physio-Chemical properties of protein using EXPASY programming.
- Multiple Sequence Alignment and BLAST on query uploaded sequence to obtain similar sequences for phylogeny
- In built protein structure visualiser through PDB code and separate locator for desired amino acid.
- Protein-Protein or Protein-Chemical interaction using STRING server API to analyse protein connections and Enrichment analysis along with GO terminologies

---

## Dependancies to install : Perform in python/conda virtual environment 

* Biopython : Biological analysis library
```sh
pip install biopython
```
* Streamlit : Python Web Development library
```sh
pip install streamlit
```
* StMol : For rendering molecules into streamlit
```sh
pip install stmol
pip install ipython_genutils
```
* Py3DMol : For rendering molecules into interactive mode python
```sh
pip install py3Dmol==2.0.0.post2
```
* Requests : To fetch API calls from developers guide 
```sh
python -m pip install requests
```
* To run the app 
```sh
streamlit run [filename]
```
---

### Links for reference 

* StMol : [StMol Github Documentation](https://github.com/napoles-uach/stmol)
* Py3DMol : [Py3DMol Github Documentation](https://www.insilicochemistry.io/tutorials/foundations/chemistry-visualization-with-py3dmol)

---