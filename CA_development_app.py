# Import Base Libraries

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Blast import NCBIWWW,NCBIXML
import streamlit as st
from stmol import showmol
from stmol import * # pip install stmol==0.0.9 , pip install ipython_genutils
import py3Dmol  # pip install py3Dmol==2.0.0.post2
import pandas as pd
import requests # # python -m pip install requests
import io  # Import the 'io' module for StringIO
from io import StringIO
from io import BytesIO


# Create a sidebar
st.sidebar.header("Sequence Operations🧬💻")
sidebar_render = st.sidebar.radio("Navigate to : ",["Home","Sequence Analysis" , "Compute Protein Structure Parameters","Sequence Alignment & BLAST","Protein Structure Visualiser","Protein-Protein and Protein-Chemical Interaction" , "Protein Modelling : ESMFold" , "Link to Our QSAR Server" , "About Us"])

# Header Main Page
if sidebar_render == "Home":
    st.title('Bioinformatics Sequence Analysis Server')
    st.write("""
    ## Welcome to the Bioinformatics Server
    ##### This server will provide base sequence analysis workflows for protein sequences in form of FASTA file format.In this project although we are considering the main workflow focus on Carbonic Anhydrases and their amino acid sequence analysis through bioinformatics analysis,this server has been build to analyse a wide spectrum of amino acid sequences in FASTA format. Some of the base features which are included into this server are :
    - FASTA Amino Acid file uploader.
    - Protein Parameters computation through ProtParam Expasy Server API.
    - BLASTP functionality and MSA of query sequence across the protein database.
    - Protein Structure Visualiser tool from PDB Code.
    - STRING server Protein-Protein and Protein-Molecule Interaction through API.(Under Development)
    - Ab initio Modelling Server through ESM Fold API.(Under Development)
    ---        
    """)

# Create Sequence Analysis Web functionality
    
if sidebar_render == "Sequence Analysis":
    st.title("Sequence Analysis Page")
    st.write("Upload your Protein sequence file in FASTA Format ⬇️")

    # Create functions for reading and decoding FASTA file and displaying it since streamlit doesnt directly render biological files like biopython CLI output

    def read_fasta_file(uploaded_file):
        fasta_file_content = uploaded_file.read().decode() # read() reads the file into streamlit and decode() decodes the bytes format into string for rendering

        return list(SeqIO.parse(StringIO(fasta_file_content),"fasta"))
    
    # Create a function to display the contents of FASTA file

    def display_fasta_file(sequence_contents):
        if len(sequence_contents) == 0:
            st.write("Sorry No Sequences Found ! 😟")
        else:
            for i,records in enumerate(sequence_contents,start=1):
                st.write(f"Sequence : {i} , ID : {records.id}")
                st.write("Your Amino Acid Sequence from the file : ")
                st.write(records.seq)
                st.write(f"Length of Sequence : {len(records.seq)}")
                st.write(f"Description of Sample : {records.description}")
                st.write("---")
        
    # Now create the file uploader 
    
    uploaded_file = st.file_uploader("Upload the FASTA file of Protein Sample 🖥️🧬" , type=["fasta"])

    if uploaded_file is not None:
        sequence_contents = read_fasta_file(uploaded_file)
        display_fasta_file(sequence_contents)
        st.success("You can copy the sequence and perform sequence operations.",icon="😍")


# Create the protein parameter analysis functionality
        
if sidebar_render == "Compute Protein Structure Parameters":
    st.title("Protein Structure Parameters Calculation Page")
    st.write("Enter your Protein sequence ⬇️")
    sequence_input = st.text_area("Enter Sequence Here. ")
    pH = st.number_input("State the pH you want to analyse your protein sequence.")
    if st.button("Compute Properties"):
        sequence_reference = ProteinAnalysis(str(sequence_input))
        st.write("<h4 style='text-align: center;'>Computed Properties of Query Protein</h4>", unsafe_allow_html=True)

        # Number of amino acids

        st.write(f"Number of each amino acids : {sequence_reference.count_amino_acids()}")

        # Compute the Molecular Weight

        st.write(f"Molecular Weight : {round(sequence_reference.molecular_weight(),2)}")

        # Compute the Aromaticity : Aromaticity of amino acids refers to the presence of an aromatic ring in their side chains, crucial for protein stability, ligand binding, and catalytic activity : The two common aromatic amino acids found in proteins are phenylalanine (Phe), tyrosine (Tyr), and tryptophan (Trp)

        st.write(f"Aromaticity : {round(sequence_reference.aromaticity(),2)}")

        # Compute the Instability index of amino acid : Value >= 40 means unstable sequence and Value < 40 means stable sequence

        stability = "The Protein is Unstable" if sequence_reference.instability_index() >= 40 else "The Protein is Stable"
        st.write(f"Instability Index: {round(sequence_reference.instability_index(), 2)} : {stability}")

        # Compute the Isoelectric point. This is the pH of molecule at net zero charge.Required for 2D electrophoresis and protein quantification

        st.write(f"Isoelectric Point : {round(sequence_reference.isoelectric_point(),2)}")

        # Compute GRAVY (Grand average of hydropathicity) : Positive values mean hydrophobic nature while Negative values mean hydrophilic nature
        charge_gravy = "Hydrophobic Nature" if sequence_reference.gravy() > 0 else "Hydrophilic Nature"
        st.write(f"GRAVY (Grand average of hydropathicity) : {round(sequence_reference.gravy(),2)} : {charge_gravy}")

        # The calculation of the molar extinction coefficient considering both cysteines (reduced) and cystines residues (Cys-Cys-bond) accounts for the absorption of light by these amino acids in a protein. Cysteine residues can exist in two forms: reduced (with a free thiol group) and oxidized (forming a disulfide bond with another cysteine residue, known as cystine). The molar extinction coefficient takes into account the absorption properties of both forms of cysteine, allowing for a more accurate estimation of the protein's light-absorbing capacity.

        st.write(f"Molar Extinction Coefficients for reduced Cysteine conformation : {sequence_reference.molar_extinction_coefficient()[0]}")
        st.write(f"Molar Extinction Coefficients for disulfid bridges Cysteine conformation : {sequence_reference.molar_extinction_coefficient()[1]}")

        # Computing the fraction of secondary structure in protein sequence.Sequentially listed as Helix , Turn , Sheet

        secondary_structure_determine = sequence_reference.secondary_structure_fraction()

        helix_fraction = secondary_structure_determine[0]
        turn_fraction = secondary_structure_determine[1]
        sheet_fraction = secondary_structure_determine[2]

        st.write(f"Percentage of Helix in amino acid sample : {round(helix_fraction*100,2)} %")
        st.write(f"Percentage of Turn in amino acid sample : {round(turn_fraction*100,2)} %")
        st.write(f"Percentage of Sheet in amino acid sample : {round(sheet_fraction*100,2)} %")

        # Making a Computing charge of protein at a given pH

        st.write(f"Charge of the query protein at pH = {pH} :  {round(sequence_reference.charge_at_pH(pH),3)}")

# Create the Sequence Alignment and BLAST functionality
        
if sidebar_render == "Sequence Alignment & BLAST":
    def run_blast(sequence_input, e_val):
        blast_result_handle = NCBIWWW.qblast("blastp", "nr", Seq(sequence_input))
        blast_records = NCBIXML.parse(blast_result_handle)
        blast_record = next(blast_records)

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_val:
                    st.write("****Alignment****")
                    st.write("▶️ Sequence:", alignment.title)
                    st.write("▶️ Length:", alignment.length)
                    # Calculate percentage identity
                    identity = (hsp.identities / hsp.align_length) * 100
                    st.write("▶️ Percentage Identity:", f"{identity:.2f}%")
                    st.write("▶️ E value:", hsp.expect)
                    st.write("▶️ Bit Score:", hsp.bits)
                    st.write("▶️ Query Sequence")
                    st.write(hsp.query[0:75] + "...")
                    st.write("▶️ Matched Sequence")
                    st.write(hsp.match[0:75] + "...")
                    st.write("▶️ Hit Sequence")
                    st.write(hsp.sbjct[0:75] + "...")
                    st.write("---")

    # Streamlit app code
    st.title("Sequence Alignment and BLAST computation page")
    st.write("Enter your Protein sequence ⬇️")
    sequence_input_blast = st.text_area("Enter Sequence Here. ")
    st.write("Enter the E Threshold Value for BLAST")
    E_VAL = st.number_input("Enter Value")

    # Submit button
    if st.button("Submit"):
        with st.spinner("Running BLAST..."):
            # Run BLAST computation
            run_blast(sequence_input_blast, E_VAL)


# Create the protein structure visualisation functionality
# Github source : https://github.com/napoles-uach/stmol           
# Py3DMol : https://www.insilicochemistry.io/tutorials/foundations/chemistry-visualization-with-py3dmol

if sidebar_render == "Protein Structure Visualiser":
    st.title("Protein Sequence Visualiser Server")
    st.write("Input the 4 character PDB code for your protein and observe the structure here! ")

    # User input features
    ## PDB Code
    PDB_Code = st.text_input("Enter the PDB Code ",value='2YLK')
    ## Display Type
    DISPLAY_TYPE = st.selectbox("Choose the display type for protein." , ["line" , "cartoon" , "cross" , "stick" , "sphere"])
    ## Residue focus
    RESIDUE_FOCUS = st.text_input("Enter a amino acid residue you want to focus on. (Eg, - ALA ). Leave Blank if not required.")
    ## Spin animation query
    SPIN_ANIMATION = st.selectbox("Do you want a spin animation on the structure" , [True , False])
    ## Background color palette
    color_pal = st.color_picker("Pick Color","#89cff0")


    def PDB_VISUALISER(PDB_Code,color_pal,DISPLAY_TYPE,SPIN_ANIMATION,RESIDUE_FOCUS):
        xyzview = py3Dmol.view(query=f'pdb:{PDB_Code}') # Creating the viewport are to display molecule

        xyzview.setStyle({DISPLAY_TYPE:{'color':'spectrum'}}) # Selecting display type along with color which molecule be displayed

        xyzview.setBackgroundColor(color_pal)  # Setting background color style

        xyzview.zoomTo() # Zoom to center

        xyzview.spin(SPIN_ANIMATION) # Spin animation

        if RESIDUE_FOCUS is not None:
            showmol(render_pdb_resn(viewer = xyzview,resn_lst = [RESIDUE_FOCUS,]),height = 500,width=800)
        else:
            showmol(xyzview, height = 500,width=800) # Render the molecule 

    if st.button("Get The Structure"):
        with st.spinner(f"Obtaining PDB Structure for {PDB_Code}"):
            PDB_VISUALISER(PDB_Code,color_pal,DISPLAY_TYPE,SPIN_ANIMATION,RESIDUE_FOCUS)


# Create Protein-Protein interaction framework and GO annotation workflow using STRING Server API
            
if sidebar_render == "Protein-Protein and Protein-Chemical Interaction":
    st.title("Protein Interaction Analysis Server")
    st.write("This server functionality provides the Network based analysis of Protein Interactions with other Proteins and Chemicals.Utilising the API from STRING Server backend the user can input their desired protein within our server and obtain the interaction network results along with their accuracy metrics and detailed GO (Gene Ontology) details and pathway enrichment analysis.")
    
    # Create input for protein name

    STRING_Protein_Query = st.text_input("Enter the accurate name of your protein eg: Rv3272")
    # Create selector for functions

    FUNCTION_QUERY = st.selectbox("What results do you want to fetch?" ,['Network Analysis Details' , 'GO Enrichment Analysis' , 'Functional Annotation of GO'])

    # Fetching The query protein interaction network

    if FUNCTION_QUERY == "Network Analysis Details":
        response = requests.get(f"https://string-db.org/api/tsv/network?identifiers={STRING_Protein_Query}")
    
    # Fetching the Enrichment GO for query and connections
        
    if FUNCTION_QUERY == "GO Enrichment Analysis":
        response = requests.get(f"https://string-db.org/api/tsv/enrichment?identifiers={STRING_Protein_Query}")

    # Fetching the annotated details for GO
        
    if FUNCTION_QUERY == "Functional Annotation of GO":
        response = requests.get(f"https://string-db.org/api/tsv/functional_annotation?identifiers={STRING_Protein_Query}")
        
    # Check if the request was successful
    if STRING_Protein_Query == "":
        st.info(f"Provide your protein")
    if response.status_code == 200:
        # Read the content as a DataFrame
        df = pd.read_csv(io.StringIO(response.text), sep='\t')

        # Display the DataFrame
        st.write(df)
    elif response.status_code == 400:
        st.error(f"Data is not available. Try Again !! Error: {response.status_code}")

    ## Fetching the network Image
    
    response = requests.get(f"https://string-db.org/api/image/network?identifiers={STRING_Protein_Query}")

    if response.status_code == 200:
        image_return =  BytesIO(response.content)
        st.header("Network Image for Protein interaction across Database")
        st.image(image_return,caption="Network")
        
    else:
        st.info(f"Provide your protein")
    
# Create the protein modelling feature through ESMFold
        
if sidebar_render == "Protein Modelling : ESMFold":
    
    def render_mol(pdb_file,SPIN_ANIMATION,RESIDUE_FOCUS):
        pdbview = py3Dmol.view()
        pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
        pdbview.addModel(pdb_file, 'pdb')
        pdbview.setBackgroundColor(color_pal)  # Setting background color style
        pdbview.zoomTo()
        pdbview.spin(SPIN_ANIMATION)
        if RESIDUE_FOCUS is not None:
            showmol(render_pdb_resn(viewer = pdbview,resn_lst = [RESIDUE_FOCUS,]),height = 500,width=800)
        else:
            showmol(pdbview, height = 500,width=800)

    st.title("Protein Structure Predictor and Modellor")
    st.header("Powered by ESM Fold API")

    AMINO_ACID_SEQUENCE_INPUT = st.text_input("Enter the amino acid sequence.Example sequence is given replace with your own!",value="AYML")

    RESIDUE_FOCUS = st.text_input("Enter a amino acid residue you want to focus on. (Eg, - ALA ). Leave Blank if not required.")
    
    ## Spin animation query
    SPIN_ANIMATION = st.selectbox("Do you want a spin animation on the structure" , [True , False])

    ## Background color palette
    color_pal = st.color_picker("Pick Color","#89cff0")

    # Disable SSL verification
    requests.packages.urllib3.disable_warnings()  # Disable script based warnings

    # Send the POST request with the query data

    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    response = requests.post(url, data=AMINO_ACID_SEQUENCE_INPUT,verify=False)

    # Check if the request was successful
    if response.status_code == 200:
        pdb_string = response.content.decode('utf-8')  # Convert the binary into text format to pass it to modeller

        st.download_button(label="Download PDB", data=pdb_string, file_name="predicted_structure.pdb", mime="chemical/x-pdb")

        st.success("Protein structure prediction successful.")
        render_mol(pdb_string,SPIN_ANIMATION,RESIDUE_FOCUS)
    else:
        st.error(f"Error: {response.status_code}, {response.reason}")

    
# Create the page to navigate to QSAR Server developed by our team
        
if sidebar_render == "Link to Our QSAR Server":
    st.title("Visit our Cheminformatics Server")
    st.write("##### Here you will be able to access the server developed by our team which specificallyd deals with In silico QSAR analysis on Carbonic Anhydrases(CA) of Mycobacterium tuberculosis bacteria.")

    st.markdown("""#### Visit the following links to access QSAR based servers : \n * QSAR Server for general analysis of Carbonic Anhydrases(CA) of Mycobacterium tuberculosis : [Carbonic Anhydrase General Analysis](https://mtb-ca-pred.streamlit.app/) \n * QSAR Server for analysis of Carbonic Anhydrases(CA) of specifically CA-1 and CA-2 of Mycobacterium tuberculosis : [Carbonic Anhydrase CA-1 and CA-2 Analysis](https://mtbca-selec-pred.streamlit.app/) """)

if sidebar_render == "About Us":
    st.title("Meet the Team of Developers and Researchers")
    st.markdown("""---""")
    st.markdown("""1. **Rajarshi Ray** : Department of Biotechnology, University of Engineering and Management Kolkata.\n --- \n 2. **Ratul Bhowmick** : Faculty of Medicine and Health Technology, Tampere University, Tampere, Finland. \n --- \n 3. **Ajay Manaithiya** : Department of Pharmaceutical Chemistry, School of Pharmaceutical Education and Research, Jamia Hamdard, New Delhi, India. \n --- \n 4. **Dr. Ashok Aspatwar** : Faculty of Medicine and Health Technology, Tampere University, Tampere, Finland. \n ---
""")