import streamlit as st
from Bio.Seq import Seq
import os
import re
from docx import Document
from PyPDF2 import PdfReader

# =======================
# UGENE-style Consensus Logic with Updated Match Points
# =======================
from pathlib import Path
import streamlit as st
from Bio.Seq import Seq

# --- Helper Functions ---
def parse_fasta_txt(txt):
    lines = txt.strip().splitlines()
    seqs = {}
    key = None
    for line in lines:
        if line.startswith(">"):
            key = line[1:].strip()
            seqs[key] = ""
        elif key:
            seqs[key] += line.strip()
    return seqs

def find_first_match_pos(s1, s3, k=5):
    for i in range(len(s1) - k + 1):
        window = s1[i:i+k]
        for j in range(len(s3) - k + 1):
            if window == s3[j:j+k]:
                return i
    return -1

def find_last_match_pos(s1, s3, threshold):
    for i in reversed(range(len(s1) - threshold)):
        if s1[i] == s3[i]:
            mismatches = sum(1 for a, b in zip(s1[i:i+threshold], s3[i:i+threshold]) if a != b)
            if mismatches >= threshold:
                return i
    return -1

def ugene_style_consensus(s1, s2, threshold=10):
    s1 = s1.strip().replace("\n", "").upper()
    s2 = s2.strip().replace("\n", "").upper()
    s2_rev = s2[::-1]
    s3 = str(Seq(s2_rev).complement())

    first_match = find_first_match_pos(s1, s3, k=5)
    last_match = find_last_match_pos(s1, s3, threshold)

    mismatch_count = -1
    a2 = ""
    if first_match != -1 and last_match != -1 and last_match > first_match:
        a2_s1 = s1[first_match:last_match+1]
        a2_s3 = s3[first_match:last_match+1]
        mismatch_count = sum(1 for a, b in zip(a2_s1, a2_s3) if a != b)
        a2 = a2_s1
    else:
        first_match = last_match = -1
        a2 = ""

    if mismatch_count > threshold:
        disclaimer = f"Number of mismatches is {mismatch_count} and it exceeds the threshold = {threshold}"
    else:
        disclaimer = None

    a1 = s1[:first_match].lower() if first_match != -1 else ""
    a3 = Seq(s3[last_match+1:]).reverse_complement().lower() if last_match != -1 else ""

    consensus = a1 + a2.upper() + str(a3)
    return consensus, s1, s2, s2_rev, s3, first_match, last_match, mismatch_count, disclaimer

# --- Streamlit App ---
st.title("IGHV Consensus Generator and Report Tool")

# SECTION 1: Consensus Generator
st.header("Section 1: UGENE-style Consensus Generator")

uploaded_txt = st.file_uploader("Upload TXT file with IGHV_F and IGHV_R reads", type=["txt"])
threshold = st.number_input("Mismatch Threshold (T)", min_value=0, max_value=100, value=10)

if uploaded_txt:
    txt_contents = uploaded_txt.read().decode("utf-8")
    seqs = parse_fasta_txt(txt_contents)
    forward = next((v for k, v in seqs.items() if k.endswith("_F")), None)
    reverse = next((v for k, v in seqs.items() if k.endswith("_R")), None)

    if forward and reverse:
        consensus, s1, s2, s2_rev, s3, first_match, last_match, mismatch_count, disclaimer = ugene_style_consensus(forward, reverse, threshold)

        st.subheader("Inputs and Intermediates")
        st.text_area("IGHV_F (Forward Read)", s1, height=100)
        st.text_area("IGHV_R (Reverse Read)", s2, height=100)
        st.text_area("Reverse of IGHV_R", s2_rev, height=100)
        st.text_area("Reverse-Complement of IGHV_R", s3, height=100)
        st.write(f"First Match Point in IGHV_F: {first_match}")
        st.write(f"Last Match Point in IGHV_F: {last_match}")
        st.write(f"Number of Mismatches: {mismatch_count}")
        if disclaimer:
            st.warning(disclaimer)

        st.subheader("Consensus Output")
        st.text_area("Consensus", consensus, height=100)

        consensus_filename = Path(uploaded_txt.name).stem + "_consensus.txt"
        st.download_button("Download Consensus", consensus, file_name=consensus_filename, mime="text/plain")
    else:
        st.error("TXT must contain both IGHV_F and IGHV_R entries.")


# =======================
# Section 2: IGHV DOCX Report Generator
# =======================
st.header("2. IGHV DOCX Report Generator")

pdf_file = st.file_uploader("Upload IgBLAST PDF", type=["pdf"])
docx_template = st.file_uploader("Upload DOCX Template", type=["docx"])
ab1_file = st.file_uploader("Upload AB1 Sequence File", type=["ab1"])

def extract_text_from_pdf(pdf_file):
    reader = PdfReader(pdf_file)
    return "\n".join(page.extract_text() for page in reader.pages if page.extract_text())

def parse_ighv_info(text):
    gene_names = re.findall(r"(IGHV\d+-\d+)", text)
    percents = re.findall(r"(\d+\.\d+)% identity", text)
    ratios = re.findall(r"\((\d+/\d+)\)", text)
    return gene_names, percents, ratios

def determine_prognosis(percent_identity, gene_name):
    mutation_rate = 100 - float(percent_identity)
    if "3-21" in gene_name:
        return "BAD"
    elif mutation_rate < 2:
        return "BAD"
    elif 2 <= mutation_rate <= 3:
        return "Borderline"
    else:
        return "GOOD"

def extract_sample_id(ab1_file):
    return os.path.splitext(ab1_file.name)[0] if ab1_file else "Sample_ID"

def replace_text_in_docx(doc, old, new):
    for p in doc.paragraphs:
        if old in p.text:
            p.text = p.text.replace(old, new)
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                if old in cell.text:
                    cell.text = cell.text.replace(old, new)

if st.button("Generate IGHV DOCX Report"):
    if pdf_file and docx_template and ab1_file:
        text = extract_text_from_pdf(pdf_file)
        gene_names, percents, ratios = parse_ighv_info(text)
        if not gene_names or not percents:
            st.error("Failed to extract IGHV info from PDF.")
        else:
            gene_name = gene_names[0]
            percent = percents[0]
            ratio = ratios[0] if ratios else ""
            prognosis = determine_prognosis(percent, gene_name)
            sample_id = extract_sample_id(ab1_file)

            doc = Document(docx_template)
            replace_text_in_docx(doc, "SAMPLE_ID", sample_id)
            replace_text_in_docx(doc, "PERCENT", f"{percent}%")
            replace_text_in_docx(doc, "PROGNOSIS", prognosis)

            for table in doc.tables:
                if table.cell(0, 0).text.strip().lower() == "gene":
                    for i in range(1, len(table.rows)):
                        table.cell(i, 0).text = ""
                        table.cell(i, 1).text = ""
                    table.cell(1, 0).text = gene_name
                    table.cell(1, 1).text = f"{percent}% ({ratio})" if ratio else f"{percent}%"

            output_docx = f"{sample_id}_IGHV_Report.docx"
            doc.save(output_docx)
            with open(output_docx, "rb") as f:
                st.download_button("Download IGHV DOCX Report", f, file_name=output_docx)
    else:
        st.error("Please upload all required files (PDF, DOCX, AB1)")
