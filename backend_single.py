import pandas as pd
import requests
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from bs4 import BeautifulSoup
import base64
import time
import streamlit as st
import random 
random_num = random.randint(0,1000)

df = pd.read_excel('Data/FPKM_Matrix(Ca).xlsx')
miRNA_df = pd.read_excel('Data/8.xlsx',header=1)
protein_df = pd.read_excel('Data/9.xlsx')
combined_data = pd.read_excel('Data/7.xlsx')
GO_df=pd.read_excel("Data/10.xlsx", sheet_name="PCGs")
cello_df=pd.read_excel("Data/13.xlsx")

df = df.reset_index(drop=True)
miRNA_df = miRNA_df.reset_index(drop=True)
protein_df = protein_df.reset_index(drop=True)
combined_data = combined_data.reset_index(drop=True)
GO_df = GO_df.reset_index(drop=True)
cello_df = cello_df.reset_index(drop=True)

def normalize_data(data):
    return data.applymap(lambda x: np.log2(x) if x > 0 else 0)

def format_sequence(seq):
    if isinstance(seq, float) and np.isnan(seq):
        return ''
    
    return '\n'.join('\t\t ' + ' '.join([seq[i:i+6] for i in range(j, min(j + 90, len(seq)), 6)]) for j in range(0, len(seq), 90))



def get_string_network_link(protein_transcript):

    api_url = "https://string-db.org/api/tsv/get_link?"

    params = {
        'identifiers': protein_transcript,
        'species': 3827,
        'format': 'json'
    }

    response = requests.get(api_url, params=params)

    if response.status_code == 200:
        return response.text.strip()
    else:
        return f"Error: {response.status_code}"

def filter_orthologs(tid):
    with open("Data/14.txt", 'r') as infile:
        lines = infile.readlines()

    filtered_data = set()

    for line in lines:
        if tid in line:
            columns = line.strip().split()
            species_a, species_b = columns[0], columns[1]
            species_pair = tuple(sorted([species_a, species_b]))
            filtered_data.add((species_pair, columns[2]))

    filtered_data_list = [(pair[0][0], pair[0][1], pair[1]) for pair in filtered_data]
    ortho_df = pd.DataFrame(filtered_data_list, columns=["Species A", "Species B", "Score"])

    return ortho_df

def filter_paralogs(tid):
    with open("Data/15.txt", 'r') as infile:
        lines = infile.readlines()

    filtered_data = set()

    for line in lines:
        if tid in line:
            columns = line.strip().split()
            species_a, species_b = columns[0], columns[1]
            species_pair = tuple(sorted([species_a, species_b]))
            filtered_data.add((species_pair, columns[2]))

    filtered_data_list = [(pair[0][0], pair[0][1], pair[1]) for pair in filtered_data]
    para_df = pd.DataFrame(filtered_data_list, columns=["Species A", "Species B", "Score"])

    return para_df

def web_driver():
    options = webdriver.ChromeOptions()
    options.add_argument("--verbose")
    options.add_argument('--no-sandbox')
    options.add_argument('--headless')
    options.add_argument('--disable-gpu')
    options.add_argument("--window-size=1920, 1200")
    options.add_argument('--disable-dev-shm-usage')
    driver = webdriver.Chrome(options=options)
    return driver

def automate_Cultivated_task(tid):
    driver = web_driver()
    driver.get("https://cegresources.icrisat.org/cicerseq/?page_id=3605")
    time.sleep(3)

    gene_id_dropdown = Select(driver.find_element(By.NAME, "select_crop"))
    gene_id_dropdown.select_by_value("cultivars")

    radio_button = driver.find_element(By.ID, "gene_snp")
    radio_button.click()

    gene_id_dropdown = Select(driver.find_element(By.NAME, "key1"))
    gene_id_dropdown.select_by_value("GeneID")

    intergenic_dropdown = Select(driver.find_element(By.NAME, "key4"))
    intergenic_dropdown.select_by_value("intergenic")

    input_field = driver.find_element(By.ID, "tmp1")
    input_field.clear()
    input_field.send_keys(tid) #Ca_00004

    search_button = driver.find_element(By.NAME, "submit")
    search_button.click()

    time.sleep(5)

    page_source = driver.page_source
    soup = BeautifulSoup(page_source, 'html.parser')
    driver.quit()

    return page_source

def automate_Wild_task(tid):
    driver = web_driver()
    driver.get("https://cegresources.icrisat.org/cicerseq/?page_id=3605")
    time.sleep(3)

    gene_id_dropdown = Select(driver.find_element(By.NAME, "select_crop"))
    gene_id_dropdown.select_by_value("wild")

    radio_button = driver.find_element(By.ID, "wgene_snp")
    radio_button.click()

    gene_id_dropdown = Select(driver.find_element(By.NAME, "key2"))
    gene_id_dropdown.select_by_value("GeneID")

    intergenic_dropdown = Select(driver.find_element(By.NAME, "key4"))
    intergenic_dropdown.select_by_value("intergenic")

    input_field = driver.find_element(By.ID, "tmp3")
    input_field.clear()
    input_field.send_keys(tid) #Ca_00004

    search_button = driver.find_element(By.NAME, "submitw")
    search_button.click()

    time.sleep(5)

    page_source = driver.page_source

    soup = BeautifulSoup(page_source, 'html.parser')
    driver.quit()

    return page_source
    
col1, col2 = st.columns(2)

def transcriptid_info(tid):
    if 'Transcript id' in df.columns and 'lncRNA' in df.columns:
        matching_row = df[df['Transcript id'] == tid]

        if not matching_row.empty:
            temp_df = df.copy()
            st.subheader("FPKM Matrix Atlas data")
            result = temp_df[temp_df['Transcript id'] == tid]
            result = result.drop(columns=['Genomic Coordinates', 'mRNA', 'lncRNA','Genomic Sequence','Transcript Sequence','Peptide Sequence','Cds Sequence','Promoter Sequence'])
            st.dataframe(result)

            st.subheader("RNA data")
            if pd.notna(matching_row['mRNA'].values[0]):
                st.write("mRNA present : Yes ( 1 )\n")
            else:
                st.write("mRNA absent : No ( 0 )\n")

            st.subheader("lncRNA data")
            if pd.notna(matching_row['lncRNA'].values[0]):
                st.write("lncRNAs present : Yes ( 1 )")
            else:
                st.write("lncRNAs absent : No ( 0 )\n")

            st.subheader("miRNA data")
            miRNA_matching_rows = miRNA_df[miRNA_df['Target_Acc.'] == tid]
            if not miRNA_matching_rows.empty:
                st.dataframe(miRNA_matching_rows)
                st.write("\n")
            else:
                st.write(f"No match found for Transcript id: {tid} in miRNA data\n")

            st.subheader("Protein and PPI data")
            protein_matching_row = protein_df[protein_df['Transcript id'] == tid]
            if not protein_matching_row.empty:
                st.dataframe(protein_matching_row)
                st.write("\n")
                protein_transcript = protein_matching_row['preferredName'].values[0]
                st.write(f"Protein Transcript for {tid}: {protein_transcript}")

                network_link = get_string_network_link(protein_transcript)
                st.write("Redirected Network URL -->", network_link)
                st.write("\n")
            else:
                st.write(f"No match found for Transcript id: {tid} in protein data\n")

            #Orthologous analysis
            st.subheader("Orthologs data")
            ortho_df = filter_orthologs(tid)
            if not ortho_df.empty:
                st.dataframe(ortho_df)
                st.write("\n")
            else:
                st.write(f"No match found for Transcript id: {tid} in Orthologs data\n")
            st.subheader("Inparalogs data")
            para_df = filter_paralogs(tid)
            if not para_df.empty:
                st.dataframe(para_df)
                st.write("\n")
            else:
                st.write(f"No match found for Transcript id: {tid} in Inparalogs data\n")
            st.write("For detailed results visit the following link -->","https://orthovenn3.bioinfotoolkits.net/result/88e9a64330ba4d64b78fc5fd9561cd64/orthologous\n")

            st.subheader("SNP Calling data")
            st.write("Result data for both Cultivated and Wild varieties will be downloaded in the form of HTML content. Click on the files to view data\n")
            try:
                # Cultivated SNP Download Button
                with col1:
                    html_Cultivated_page_source = automate_Cultivated_task(tid)
                    b64_html = base64.b64encode(html_Cultivated_page_source.encode()).decode()  # Convert to base64
                    html_href = f'<a href="data:text/html;base64,{b64_html}" download="{tid}_Cultivated_SNP.html">Download Cultivated SNP as .html</a>'
                    st.markdown(html_href, unsafe_allow_html=True)
                # Wild SNP Download Button
                with col2:
                    html_wild_page_source = automate_Wild_task(tid)
                    b64_html2 = base64.b64encode(html_wild_page_source.encode()).decode()  # Convert to base64
                    html_href2 = f'<a href="data:text/html;base64,{b64_html2}" download="{tid}_Wild_SNP.html">Download Wild SNP as .html</a>'
                    st.markdown(html_href2, unsafe_allow_html=True)
            except Exception as e:
                st.write("Error ! Error ! Error !")
                st.write("Unable to fetch data from the server. Please try again later -->","https://cegresources.icrisat.org/cicerseq/?page_id=3605\n")
            st.subheader("GO and KEGG data")
            GO_matching_row = GO_df[GO_df['Transcript id'] == tid]
            if not GO_matching_row.empty:
                st.dataframe(GO_matching_row)
                st.write("\n")
            else:
                st.write(f"No match found for Transcript id: {tid} in GO KEGG data\n")

            st.subheader("Cellular Localisation data")
            cello_matching_row = cello_df[cello_df['Transcript id'] == tid]
            if not cello_matching_row.empty:
                cello_matching_row = cello_matching_row.head(1)
                cello_matching_row = cello_matching_row.drop(columns=["#Combined:"])
                st.dataframe(cello_matching_row)
                st.write("\n")
            else:
                st.write(f"No match found for Transcript id: {tid} in Cellular Protein Localisation data\n")

            st.subheader("Sequences data")
            cds_code = format_sequence(matching_row['Cds Sequence'].values[0])
            peptide_code = format_sequence(matching_row['Peptide Sequence'].values[0])
            transcript_code = format_sequence(matching_row['Transcript Sequence'].values[0])
            gene_code = format_sequence(matching_row['Genomic Sequence'].values[0])
            promote_code = format_sequence(matching_row['Promoter Sequence'].values[0])

            # Display as code block with copy functionality
            st.markdown("### CDS Sequence")
            st.code(cds_code, language="text")
            st.markdown("### Peptide Sequence")
            st.code(peptide_code, language="text")
            st.markdown("### Transcript Sequence")
            st.code(transcript_code, language="text")
            st.markdown("### Genomic Sequence")
            st.code(gene_code, language="text")
            st.markdown("### Promoter Sequence")
            st.code(promote_code, language="text")

            header = f">{tid}|{tid}"
            promote_file = f"{header}\n{promote_code}\n"
            b64 = base64.b64encode(promote_file.encode()).decode()  # Convert to base64
            href = f'<a href="data:text/plain;base64,{b64}" download="{tid}_promoter_sequence.txt">Download Promoter Sequence as .txt</a>'
            st.markdown(href, unsafe_allow_html=True)
            #st.download_button( label="Download Promoter Sequence as .txt", data=promote_file, file_name=f"{tid}_promoter_sequence.txt", mime="text/plain" )
            st.write("Paste the promoter sequence on the following link to get promoter region analysis!")
            st.write("https://bioinformatics.psb.ugent.be/webtools/plantcare/html/search_CARE_onCluster.html\n")

        else:
            st.write("Transcript ID not found\n")
    else:
        st.write("...Error...\n")

def user_input_menu(tid):
        transcriptid_info(tid)
        if tid in combined_data['Transcript id'].values:
            st.subheader("Model Prediction")
            resultant_value = combined_data[combined_data['Transcript id'] == tid]['Resultant'].values[0]
            st.write(f"Stage/Tissue Group Concerned {tid}: {resultant_value}\n")
        else:
            st.subheader("Model Prediction")
            st.write("Expression Status : normal  ( no particular tissue/stage favoured ) 0 \n")