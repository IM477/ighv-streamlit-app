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
from Bio.Seq import Seq
from Bio import pairwise2
from io import StringIO

# Reverse complement function
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

# Function to parse FASTA-style input
def parse_fasta(txt):
    records = {}
    key = None
    seq = []
    for line in txt.strip().splitlines():
        line = line.strip()
        if line.startswith('>'):
            if key and seq:
                records[key] = ''.join(seq)
            key = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    if key and seq:
        records[key] = ''.join(seq)
    return records

# Function to calculate alignment and match points
def local_alignment_match_points(s1, s3, mismatch_threshold=10):
    alignments = pairwise2.align.localms(s1, s3, 2, -1, -0.5, -0.1)
    if not alignments:
        return -1, -1, -1, "", "", True

    best = alignments[0]
    aligned_s1 = best.seqA
    aligned_s3 = best.seqB
    first = best.start
    last = best.end - 1
    mismatches = sum(1 for a, b in zip(aligned_s1, aligned_s3) if a != b)

    alignment_str = ""
    for a, b in zip(aligned_s1, aligned_s3):
        if a == b:
            alignment_str += "|"
        else:
            alignment_str += " "

    return first, last, mismatches, aligned_s1, aligned_s3, mismatches > mismatch_threshold

# Streamlit app
st.title("UGENE-Style Consensus Viewer with Alignment")

uploaded_file = st.file_uploader("Upload TXT file (FASTA format)", type="txt")

if uploaded_file:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    fasta_data = stringio.read()
    sequences = parse_fasta(fasta_data)

    forward_label = next((k for k in sequences if k.endswith("_F")), None)
    reverse_label = next((k for k in sequences if k.endswith("_R")), None)

    if forward_label and reverse_label:
        s1 = sequences[forward_label].upper()
        s2 = sequences[reverse_label].upper()
        s3 = reverse_complement(s2)

        first, last, mismatches, aligned_s1, aligned_s3, exceeds_threshold = local_alignment_match_points(s1, s3)

        st.subheader("Input Sequences")
        st.text_area("Forward Read (s1)", s1, height=150)
        st.text_area("Reverse Read (s2)", s2, height=150)
        st.text_area("Reverse Complement (s3)", s3, height=150)

        st.subheader("Local Alignment")
        if first != -1:
            st.text("s1:\n" + aligned_s1)
            st.text("   \n" + ''.join(['|' if a == b else ' ' for a, b in zip(aligned_s1, aligned_s3)]))
            st.text("s3:\n" + aligned_s3)
        else:
            st.warning("No alignment found between forward and reverse-complement.")

        st.subheader("Match Info")
        st.write(f"**First Match Point in IGHV_F:** {first}")
        st.write(f"**Last Match Point in IGHV_F:** {last}")
        st.write(f"**Number of Mismatches:** {mismatches}")

        if exceeds_threshold:
            st.error("⚠️ Mismatch count exceeds threshold – review required.")
        else:
            st.success("✅ Mismatch count is within threshold.")
    else:
        st.error("Could not find both _F and _R labeled sequences in the file.")

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
