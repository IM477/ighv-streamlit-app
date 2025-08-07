import streamlit as st
from Bio.Seq import Seq
import os
import re
from docx import Document
from PyPDF2 import PdfReader

# =======================
# UGENE-style Consensus Logic with Updated Match Points
# =======================
import streamlit as st
from pathlib import Path
from Bio.Seq import Seq
from io import StringIO
import base64

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_first_last_match(s1, s3, match_length=6, mismatch_threshold=10):
    first_match_index = -1
    last_match_index = -1

    # Find first match: left to right
    for i in range(len(s1) - match_length + 1):
        window_s1 = s1[i:i + match_length]
        for j in range(len(s3) - match_length + 1):
            if window_s1 == s3[j:j + match_length]:
                first_match_index = i
                break
        if first_match_index != -1:
            break

    # Find last match: right to left
    for i in range(len(s1) - match_length, -1, -1):
        window_s1 = s1[i:i + match_length]
        for j in range(len(s3) - match_length, -1, -1):
            if window_s1 == s3[j:j + match_length]:
                last_match_index = i + match_length - 1
                break
        if last_match_index != -1:
            break

    mismatches = -1
    exceeds_threshold = False
    if first_match_index != -1 and last_match_index != -1 and first_match_index < last_match_index:
        region_s1 = s1[first_match_index:last_match_index + 1]
        region_s3 = s3[first_match_index:last_match_index + 1]
        mismatches = sum(1 for a, b in zip(region_s1, region_s3) if a != b)
        if mismatches > mismatch_threshold:
            exceeds_threshold = True

    return first_match_index, last_match_index, mismatches, exceeds_threshold

def compute_consensus(s1, s3, first, last):
    if first == -1 or last == -1 or first >= last:
        return ""
    consensus_prefix = s1[:last + 1]
    suffix_s3 = s3[last + 1:]
    consensus_suffix = reverse_complement(suffix_s3).lower()
    return consensus_prefix + consensus_suffix

def parse_txt_fasta(content):
    seqs = {}
    lines = content.strip().splitlines()
    current_key = ""
    for line in lines:
        if line.startswith(">"):
            current_key = line[1:].strip()
            seqs[current_key] = ""
        else:
            seqs[current_key] += line.strip()
    return seqs

def download_link(content, filename, label):
    b64 = base64.b64encode(content.encode()).decode()
    href = f'<a href="data:file/txt;base64,{b64}" download="{filename}">{label}</a>'
    return href

# Streamlit App

st.title("IGHV Consensus Generator (UGENE-style Reverse-Complement Matching)")

st.markdown("### Section 1: Reverse Complement Matching with Threshold")

uploaded_file = st.file_uploader("Upload IGHV .txt file with forward and reverse reads", type=["txt"], key="txt_upload")

T = st.number_input("Mismatch Threshold (T)", min_value=0, max_value=100, value=10)

if uploaded_file:
    content = uploaded_file.read().decode()
    seqs = parse_txt_fasta(content)

    forward = next((v for k, v in seqs.items() if k.endswith("_F")), None)
    reverse = next((v for k, v in seqs.items() if k.endswith("_R")), None)

    if forward and reverse:
        s1 = forward.strip().replace("\n", "").upper()
        s2 = reverse.strip().replace("\n", "").upper()
        s2_rev = s2[::-1]
        s3 = reverse_complement(s2)

        first, last, mismatches, exceeds = find_first_last_match(s1, s3, match_length=6, mismatch_threshold=T)
        consensus = compute_consensus(s1, s3, first, last)

        st.subheader("Inputs and Intermediate Outputs")
        st.text_area("IGHV_F (s1)", s1, height=100)
        st.text_area("IGHV_R (s2)", s2, height=100)
        st.text_area("Reverse of IGHV_R (s2_rev)", s2_rev, height=100)
        st.text_area("Reverse-Complement of IGHV_R (s3)", s3, height=100)

        st.write(f"**First Match Point in IGHV_F:** {first}")
        st.write(f"**Last Match Point in IGHV_F:** {last}")
        st.write(f"**Number of Mismatches:** {mismatches}")

        if exceeds:
            st.warning(f"Number of mismatches is {mismatches}, which exceeds the threshold = {T}")
        elif first == -1 or last == -1:
            st.error("No valid matching region found between s1 and reverse-complement of s2.")
        else:
            st.success("Consensus generated successfully.")

            st.subheader("Visual Alignment (UGENE-style)")

            align_display = ""
            for i in range(first, last + 1):
                base1 = s1[i]
                base2 = s3[i]
                marker = "|" if base1 == base2 else "."
                align_display += f"{base1} {marker} {base2}\n"
            st.text(align_display)

            st.subheader("Consensus Output")
            st.text_area("Final Consensus", consensus, height=100)

            filename = uploaded_file.name.replace(".txt", "_consensus.txt")
            st.markdown(download_link(consensus, filename, "Download Consensus Output"), unsafe_allow_html=True)

    else:
        st.error("Could not find both forward (_F) and reverse (_R) sequences in uploaded file.")

# Section 2: IGHV DOCX Report Generator (Unchanged from your existing code)
st.markdown("---")
st.markdown("### Section 2: IGHV DOCX Report Generator (Unchanged)")
st.info("This section remains unchanged. Paste your existing logic below this line.")


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
