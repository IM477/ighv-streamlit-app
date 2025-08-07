import streamlit as st
from Bio.Seq import Seq
import os
import re
from docx import Document
from PyPDF2 import PdfReader

# =======================
# UGENE-style Consensus Logic with Updated Match Points
# =======================
def ugene_style_consensus(forward, reverse, threshold=10):
    s1 = forward.strip().replace("\n", "").upper()
    s2 = reverse.strip().replace("\n", "").upper()
    s2_rev = s2[::-1]
    s3 = str(Seq(s2_rev).complement())

    # Find first match point: s1[i:i+5] == s3[i:i+5]
    first_match_pos = -1
    for i in range(len(s1) - 4):
        if s1[i:i+5] == s3[i:i+5]:
            first_match_pos = i
            break

    # Find last match point: s1[j:j+T] != s3[j:j+T]
    last_match_pos = -1
    for j in reversed(range(len(s1) - threshold)):
        if s1[j:j+threshold] != s3[j:j+threshold]:
            last_match_pos = j
            break

    if first_match_pos == -1 or last_match_pos == -1 or first_match_pos >= last_match_pos:
        return None, s1, s2, s2_rev, s3, first_match_pos, last_match_pos, None, "No valid match found"

    matched_s1 = s1[first_match_pos:last_match_pos + 1]
    matched_s3 = s3[first_match_pos:last_match_pos + 1]

    mismatches = sum(1 for a, b in zip(matched_s1, matched_s3) if a != b)

    disclaimer = None
    if mismatches > threshold:
        disclaimer = f"Number of mismatches = {mismatches}, which exceeds threshold = {threshold}"

    a1 = s1[:last_match_pos + 1]
    suffix_s3 = s3[last_match_pos + 1:]
    a3 = str(Seq(suffix_s3).reverse_complement())

    consensus = a1 + a3

    return consensus, s1, s2, s2_rev, s3, first_match_pos, last_match_pos, mismatches, disclaimer


# =======================
# Streamlit App UI
# =======================
st.set_page_config(layout="wide")
st.title("IGHV UGENE-style Consensus Generator")

# Section 1: Consensus Generator
st.header("1. Generate Consensus from Forward/Reverse TXT")
uploaded_txt = st.file_uploader("Upload Forward + Reverse Reads TXT File", type=["txt"], key="txt_file")
threshold = st.number_input("Mismatch Threshold (T)", min_value=1, value=10, step=1)

if uploaded_txt:
    contents = uploaded_txt.read().decode("utf-8").splitlines()
    seqs = {}
    current_key = None
    for line in contents:
        if line.startswith(">"):
            current_key = line[1:].strip()
            seqs[current_key] = ""
        elif current_key:
            seqs[current_key] += line.strip()

    forward = next((v for k, v in seqs.items() if k.endswith("_F")), None)
    reverse = next((v for k, v in seqs.items() if k.endswith("_R")), None)

    if forward and reverse:
        result = ugene_style_consensus(forward, reverse, threshold)
        if result is None:
            st.error("No valid match region found. Cannot generate consensus.")
        else:
            consensus, s1, s2, s2_rev, s3, first_pos, last_pos, mismatches, disclaimer = result

            st.subheader("Inputs and Intermediates")
            st.text_area("Forward Read (s1)", s1, height=100)
            st.text_area("Reverse Read (s2)", s2, height=100)
            st.text_area("Reverse of s2 (s2_rev)", s2_rev, height=100)
            st.text_area("Reverse Complement of s2 (s3)", s3, height=100)

            st.subheader("Matching Info")
            st.markdown(f"**First Match Position in s1:** {first_pos}")
            st.markdown(f"**Last Match Position in s1:** {last_pos}")
            st.markdown(f"**Mismatches between matched region of s1 and s3:** {mismatches}")
            if disclaimer:
                st.warning(disclaimer)

            st.subheader("Final Consensus")
            st.text_area("Consensus Sequence", consensus, height=100)

            output_filename = os.path.splitext(uploaded_txt.name)[0] + "_consensus.txt"
            st.download_button("Download Consensus", consensus, file_name=output_filename, mime="text/plain")
    else:
        st.error("Missing forward (_F) or reverse (_R) read in TXT file.")


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
