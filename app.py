import streamlit as st
from io import BytesIO
from docx import Document
import PyPDF2
import re
from datetime import datetime
from Bio import pairwise2

# ------------------------
# UGENE-style Consensus Logic
# ------------------------

def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def ugene_style_consensus(forward, reverse):
    s1 = forward.strip().replace("\n", "").upper()
    s2 = reverse.strip().replace("\n", "").upper()
    s2_rev = s2[::-1]
    s3 = reverse_complement(s2)

    alignments = pairwise2.align.localms(s1, s3, 2, -1, -5, -0.5)
    if not alignments:
        return None, s1, s2, s2_rev, s3, ""

    best = alignments[0]
    aligned_s1, aligned_s3, score, start, end = best

    consensus = []
    for b1, b3 in zip(aligned_s1, aligned_s3):
        if b1 == b3 and b1 != "-":
            consensus.append(b1.upper())
        elif b3 == "-":
            consensus.append(b1.lower())
        elif b1 == "-":
            consensus.append(b3.lower())
        else:
            consensus.append(b1.lower())  # mismatch

    a1 = s1[:start].lower()
    a2 = "".join(consensus)
    a3_raw = s3[end:]
    a3 = reverse_complement(a3_raw).lower()  # convert to s2 orientation

    if not a2.strip():
        return "NO_MATCH", s1, s2, s2_rev, s3, ""

    s4 = a1 + a2 + a3
    return s4, s1, s2, s2_rev, s3, a2

def parse_fasta(text):
    seqs = {}
    current_label = None
    for line in text.strip().splitlines():
        line = line.strip()
        if line.startswith(">"):
            current_label = line[1:]
            seqs[current_label] = ""
        elif current_label:
            seqs[current_label] += line.upper()
    return seqs

# ------------------------
# IGHV PDF extraction
# ------------------------
def extract_from_pdf(pdf_bytes):
    pdf_reader = PyPDF2.PdfReader(BytesIO(pdf_bytes))
    text = "".join(page.extract_text() or "" for page in pdf_reader.pages)

    gene_names_list = []
    percent_identity_str = None
    ratio_str = None

    pattern_section = re.compile(r"Sequences\s+pr.*?alignmen", re.IGNORECASE | re.DOTALL)
    section_match = pattern_section.search(text)
    if section_match:
        start_pos = section_match.end()
        lines_after_section = text[start_pos:].strip().splitlines()
        for line in lines_after_section:
            for word in line.strip().split():
                if word.startswith("IGHV"):
                    gene_names_list.append(word)
            if gene_names_list:
                break

    pattern_total = r"Total\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+([\d\.]+)"
    match_total = re.search(pattern_total, text)
    if match_total:
        length = int(match_total.group(1))
        matches = int(match_total.group(2))
        identity = float(match_total.group(3))
        percent_identity_str = f"{identity}%"
        ratio_str = f"{matches}/{length}"

    return gene_names_list, percent_identity_str, ratio_str

# ------------------------
# DOCX Report Generator
# ------------------------
def generate_docx_report(
    gene_names_list,
    percent_identity_str,
    ratio_str,
    template_bytes,
    sample_id_text
):
    gene_name_str = ", ".join(str(g).strip() for g in gene_names_list)
    percent_identity = float(percent_identity_str.strip('%'))
    mutation_percent = round(100 - percent_identity, 1)
    mutation_percent_str = f"{mutation_percent}%"

    if any("3-21" in g for g in gene_names_list):
        prognosis_text = (
            "BAD\n"
            "Note: CLL expressing the IGHV 3-21 variable region gene segment "
            "have a poorer prognosis regardless of IGHV mutation status."
        )
        prognosis_single_line = prognosis_text
    else:
        if mutation_percent < 2.0:
            prognosis_text = "BAD"
        elif 2.1 <= mutation_percent <= 3.0:
            prognosis_text = "Borderline with intermediate clinical course"
        else:
            prognosis_text = "GOOD"
        prognosis_single_line = prognosis_text

    doc = Document(BytesIO(template_bytes))

    for section in doc.sections:
        for para in section.header.paragraphs:
            if "Sample ID" in para.text:
                para.text = f"Sample ID: {sample_id_text}"

    for para in doc.paragraphs:
        if "Sample ID" in para.text:
            para.text = f"Sample ID: {sample_id_text}"
        elif "Mutation detected:" in para.text:
            para.text = f"Mutation detected: {mutation_percent_str}"
        elif "Clinical Prognosis" in para.text:
            idx = doc.paragraphs.index(para)
            if idx + 1 < len(doc.paragraphs):
                doc.paragraphs[idx + 1].text = prognosis_text

    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                if "Sample ID" in cell.text:
                    cell.text = f"Sample ID: {sample_id_text}"
                elif "Mutation detected:" in cell.text:
                    row.cells[1].text = mutation_percent_str
                elif "Clinical Prognosis:" in cell.text:
                    row.cells[1].text = prognosis_single_line

    for table in doc.tables:
        first_row = [cell.text.strip() for cell in table.rows[0].cells]
        if "Name" in first_row and "Percentage Identity Detected" in first_row:
            for i in range(len(table.rows) - 1, 0, -1):
                table._tbl.remove(table.rows[i]._tr)
            row_cells = table.add_row().cells
            row_cells[0].text = gene_name_str
            row_cells[1].text = percent_identity_str
            row_cells[2].text = mutation_percent_str
            row_cells[3].text = ratio_str
            break

    output = BytesIO()
    doc.save(output)
    output.seek(0)
    return output

# ------------------------
# STREAMLIT APP
# ------------------------
st.set_page_config(page_title="IGHV Report Generator", layout="centered")
st.title("IGHV Report Generator")
st.caption("*UGENE-style consensus generator + IGHV DOCX reporter*")

# ------------------------
# SECTION 1: Consensus Generator
# ------------------------
st.header("1. UGENE-style Consensus Generator")
cap3_input_file = st.file_uploader("Upload FASTA File with Forward (_F) and Reverse (_R) Reads", type=["txt", "fasta"])

if cap3_input_file:
    content = cap3_input_file.read().decode("utf-8")
    seqs = parse_fasta(content)
    forward = next((v for k, v in seqs.items() if k.endswith("_F")), None)
    reverse = next((v for k, v in seqs.items() if k.endswith("_R")), None)

    if forward and reverse:
        result, s1, s2, s2_rev, s3, match = ugene_style_consensus(forward, reverse)

        st.subheader("Inputs & Intermediate Steps")
        st.code(f"Forward Read (s1):\n{s1}", language="text")
        st.code(f"Reverse Read (s2):\n{s2}", language="text")
        st.code(f"Reverse of s2:\n{s2_rev}", language="text")
        st.code(f"Reverse Complement of s2 (s3):\n{s3}", language="text")

        if result == "NO_MATCH":
            st.error("⚠️ No matching region found. Manual intervention needed.")
        else:
            st.code(f"Matching Region (uppercase in output):\n{match}", language="text")
            st.success("UGENE-style Consensus Generated")
            st.code(result, language="text")

            st.download_button(
                label="Download Consensus (.txt)",
                data=BytesIO(result.encode("utf-8")),
                file_name="consensus_output.txt",
                mime="text/plain"
            )
    else:
        st.error("Could not find both _F and _R sequences in file.")

# ------------------------
# SECTION 2: IGHV Report Generator
# ------------------------
st.header("2. IGHV DOCX Report Generator")

pdf_file = st.file_uploader("PDF file (IgBLAST results)", type="pdf")
docx_template_file = st.file_uploader("Word Template (.docx)", type="docx")
ab1_file = st.file_uploader(".ab1 Sequence File", type="ab1")

if st.button("Generate IGHV Report"):
    if pdf_file and docx_template_file and ab1_file:
        ab1_filename = ab1_file.name
        sample_id = ab1_filename.split("(")[0].strip()
        date_str = datetime.now().strftime("%d-%m-%Y")
        final_name = f"{sample_id} ({date_str})"
        output_filename = f"{final_name}.docx"

        st.info(f"Sample ID: {sample_id}")
        st.info(f"Output file name: {output_filename}")

        pdf_bytes = pdf_file.read()
        gene_names_list, percent_identity_str, ratio_str = extract_from_pdf(pdf_bytes)

        st.write("Gene Names:", gene_names_list)
        st.write("Percent Identity:", percent_identity_str)
        st.write("Ratio:", ratio_str)

        if gene_names_list and percent_identity_str and ratio_str:
            docx_bytes = generate_docx_report(
                gene_names_list=gene_names_list,
                percent_identity_str=percent_identity_str,
                ratio_str=ratio_str,
                template_bytes=docx_template_file.read(),
                sample_id_text=final_name
            )

            st.download_button(
                label="Download Report",
                data=docx_bytes,
                file_name=output_filename,
                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document"
            )
        else:
            st.error("Could not extract all required fields from PDF.")
    else:
        st.warning("Please upload all three files before generating the report.")
