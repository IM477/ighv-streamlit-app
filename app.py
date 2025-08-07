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
from Bio import pairwise2

# ---------- Helper Functions ----------

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def parse_fasta(fasta_text: str) -> dict:
    lines = fasta_text.strip().splitlines()
    seqs = {}
    header = None
    for line in lines:
        if line.startswith(">"):
            header = line[1:].strip()
            seqs[header] = ""
        else:
            seqs[header] += line.strip()
    return seqs

def ugene_style_consensus(forward: str, reverse: str, threshold: int):
    s1 = forward.strip().replace("\n", "").upper()
    s2 = reverse.strip().replace("\n", "").upper()
    s2_rev = s2[::-1]
    s3 = reverse_complement(s2)

    # Local alignment
    alignments = pairwise2.align.localms(s1, s3, 2, -1, -5, -0.5)
    if not alignments:
        return "", s1, s2, s2_rev, s3, "", -1, -1, -1, "No alignment found.", "", "", ""

    best_alignment = alignments[0]
    aligned_s1 = best_alignment.seqA
    aligned_s3 = best_alignment.seqB
    start = best_alignment.start
    end = best_alignment.end

    aligned_s1_region = aligned_s1[start:end]
    aligned_s3_region = aligned_s3[start:end]
    mismatch_count = sum(1 for a, b in zip(aligned_s1_region, aligned_s3_region) if a != b)
    match_line = ''.join('|' if a == b else ' ' for a, b in zip(aligned_s1, aligned_s3))

    # Find first match point: 6-base exact match
    first_match = -1
    for i in range(len(s1) - 5):
        if s1[i:i+6] == s3[i:i+6]:
            first_match = i
            break

    # Find last match point: where threshold-length window has mismatches ≥ threshold
    last_match = -1
    for i in range(len(s1) - threshold):
        window_s1 = s1[i:i+threshold]
        window_s3 = s3[i:i+threshold]
        mismatches = sum(1 for a, b in zip(window_s1, window_s3) if a != b)
        if mismatches >= threshold:
            last_match = i + threshold - 1
            break

    a2 = s1[first_match:last_match+1] if first_match != -1 and last_match != -1 else ""
    a1 = s1[:first_match].lower() if first_match != -1 else ""
    a3 = s3[last_match+1:].lower() if last_match != -1 else ""
    consensus = a1 + a2.upper() + a3

    disclaimer = ""
    if mismatch_count > threshold:
        disclaimer = f"⚠️ Mismatches ({mismatch_count}) exceed threshold ({threshold})"

    return consensus, s1, s2, s2_rev, s3, a2, first_match, last_match, mismatch_count, disclaimer, aligned_s1, match_line, aligned_s3

# ---------- Streamlit App ----------

st.set_page_config(page_title="IGHV Consensus Generator", layout="wide")
st.title("🧬 IGHV UGENE-style Consensus Generator")

st.markdown("### 🔹 Section 1: Generate Consensus")

uploaded_file = st.file_uploader("Upload a TXT file (FASTA format with _F and _R reads)", type=["txt"])
threshold = st.number_input("Mismatch Threshold (T)", min_value=1, max_value=50, value=10)

if uploaded_file:
    content = uploaded_file.read().decode("utf-8")
    seqs = parse_fasta(content)
    forward = next((v for k, v in seqs.items() if k.endswith("_F")), None)
    reverse = next((v for k, v in seqs.items() if k.endswith("_R")), None)

    if forward and reverse:
        consensus, s1, s2, s2_rev, s3, a2, first_match, last_match, mismatch_count, disclaimer, aligned_s1, match_line, aligned_s3 = ugene_style_consensus(forward, reverse, threshold)

        st.subheader("🔍 Intermediate Results")
        st.text_area("Forward Read (s1)", s1, height=100)
        st.text_area("Reverse Read (s2)", s2, height=100)
        st.text_area("Reverse of s2", s2_rev, height=100)
        st.text_area("Reverse-Complement of s2 (s3)", s3, height=100)
        st.text_area("Matched Region (a2)", a2, height=60)

        st.write(f"🧭 First Match Point in IGHV_F: {first_match}")
        st.write(f"🔚 Last Match Point in IGHV_F: {last_match}")
        st.write(f"❌ Number of Mismatches in Alignment: {mismatch_count}")
        if disclaimer:
            st.warning(disclaimer)

        st.subheader("🧩 Visual Alignment")
        st.text("Aligned Forward Read (s1)")
        st.code(aligned_s1)
        st.text("Match Line")
        st.code(match_line)
        st.text("Aligned Reverse Complement (s3)")
        st.code(aligned_s3)

        st.subheader("✅ Final Consensus")
        st.code(consensus, language="text")

        # Download button
        output_file = Path(uploaded_file.name).with_suffix(".consensus.txt")
        st.download_button("📥 Download Consensus", consensus, file_name=output_file.name, mime="text/plain")

st.markdown("---")
st.markdown("### 🔹 Section 2: IGHV DOCX Report Generator")
st.info("This section is unchanged and handles IGHV reporting.")

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
