import streamlit as st
from io import BytesIO
from docx import Document
import PyPDF2
import re
from datetime import datetime
import subprocess
import tempfile
from pathlib import Path
import os

# ============================================================
# CAP3 Function
# ============================================================
def run_cap3_and_get_consensus(raw_txt):
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / "reads.txt"
        with open(input_path, "wb") as f:
            f.write(raw_txt)

        subprocess.run(
            ["cap3", str(input_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        ace_path = str(input_path) + ".cap.ace"

        if not os.path.exists(ace_path):
            return None, "CAP3 failed: No .ace file produced"

        consensus = ""
        in_sequence = False
        with open(ace_path, 'r') as f:
            for line in f:
                if line.startswith("CO "):
                    in_sequence = True
                    continue
                if in_sequence:
                    if line.startswith("BQ"):
                        break
                    consensus += line.strip()

        return (consensus.strip(), None) if consensus else (None, "No consensus found in .ace")

# ============================================================
# IGHV PDF Parser
# ============================================================
def extract_from_pdf(pdf_bytes):
    pdf_reader = PyPDF2.PdfReader(BytesIO(pdf_bytes))
    text = ""
    for page in pdf_reader.pages:
        text += "\n" + page.extract_text()

    gene_names_list = []
    percent_identity_str = None
    ratio_str = None

    pattern_section = re.compile(r"Sequences\s+pr.*?alignmen", re.IGNORECASE | re.DOTALL)
    section_match = pattern_section.search(text)
    if section_match:
        start_pos = section_match.end()
        remaining_text = text[start_pos:]
        lines_after_section = remaining_text.strip().splitlines()
        for line in lines_after_section:
            words = line.strip().split()
            for word in words:
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

# ============================================================
# IGHV DOCX Report Generator
# ============================================================
def generate_docx_report(
    gene_names_list,
    percent_identity_str,
    ratio_str,
    template_bytes,
    sample_id_text,
    consensus_sequence=None
):
    gene_name_str = ", ".join(gene_names_list)
    percent_identity = float(percent_identity_str.strip('%'))
    mutation_percent = round(100 - percent_identity, 1)
    mutation_percent_str = f"{mutation_percent}%"

    if any("3-21" in g for g in gene_names_list):
        prognosis_text = "BAD\nNote: IGHV 3-21 has poor prognosis regardless of mutation."
        prognosis_single_line = "BAD"
    elif mutation_percent < 2.0:
        prognosis_text = prognosis_single_line = "BAD"
    elif 2.1 <= mutation_percent <= 3.0:
        prognosis_text = prognosis_single_line = "Borderline with intermediate clinical course"
    else:
        prognosis_text = prognosis_single_line = "GOOD"

    doc = Document(BytesIO(template_bytes))

    for section in doc.sections:
        for para in section.header.paragraphs:
            if "Sample ID" in para.text:
                para.text = f"Sample ID: {sample_id_text}"

    for para in doc.paragraphs:
        if "Sample ID" in para.text:
            para.text = f"Sample ID: {sample_id_text}"
        if "Mutation detected:" in para.text:
            para.text = f"Mutation detected: {mutation_percent_str}"
        if "Clinical Prognosis" in para.text:
            idx = doc.paragraphs.index(para)
            if idx + 1 < len(doc.paragraphs):
                doc.paragraphs[idx + 1].text = prognosis_text

    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                if "Sample ID" in cell.text:
                    cell.text = f"Sample ID: {sample_id_text}"
                if "Mutation detected:" in cell.text:
                    row.cells[1].text = mutation_percent_str
                if "Clinical Prognosis:" in cell.text:
                    row.cells[1].text = prognosis_single_line

        headers = [cell.text.strip() for cell in table.rows[0].cells]
        if "Name" in headers and "Percentage Identity Detected" in headers:
            for i in range(len(table.rows) - 1, 0, -1):
                table._tbl.remove(table.rows[i]._tr)
            row_cells = table.add_row().cells
            row_cells[0].text = gene_name_str
            row_cells[1].text = percent_identity_str
            row_cells[2].text = mutation_percent_str
            row_cells[3].text = ratio_str

    if consensus_sequence:
        doc.add_paragraph("Consensus Sequence (CAP3):")
        doc.add_paragraph(consensus_sequence)

    output = BytesIO()
    doc.save(output)
    output.seek(0)
    return output

# ============================================================
# STREAMLIT APP
# ============================================================
st.set_page_config(page_title="IGHV Report Generator", layout="centered")
st.title("IGHV Report Generator")
st.caption("*A Wobble Base Bioresearch proprietary software*")

# ------------------ SECTION 1: CAP3 STANDALONE ------------------
st.header("🧬 CAP3 Consensus Generator (UGENE-style)")
cap3_input_file = st.file_uploader("Upload FASTA (.txt) with 2 reads", type=["txt"], key="cap3")

if cap3_input_file:
    if st.button("Generate Consensus from CAP3"):
        consensus, error = run_cap3_and_get_consensus(cap3_input_file.read())
        if consensus:
            st.success("Consensus generated successfully.")
            st.code(consensus, language="text")

            # Prepare downloadable .txt file
            cap3_output = BytesIO()
            cap3_output.write(consensus.encode())
            cap3_output.seek(0)

            st.download_button(
                label="Download Consensus (.txt)",
                data=cap3_output,
                file_name="cap3_consensus.txt",
                mime="text/plain"
            )
        else:
            st.error(f"CAP3 failed: {error}")

st.markdown("---")

# ------------------ SECTION 2: IGHV REPORT GENERATOR ------------------
st.header("🧾 IGHV Mutation Report Generator")

pdf_file = st.file_uploader("Upload IgBLAST PDF", type="pdf")
docx_template_file = st.file_uploader("Upload Word Template (.docx)", type="docx")
ab1_file = st.file_uploader("Upload .ab1 Sequence File", type="ab1")

if st.button("Generate Report"):
    if pdf_file and docx_template_file and ab1_file:
        sample_id = ab1_file.name.split("(")[0].strip()
        date_str = datetime.now().strftime("%d-%m-%Y")
        final_name = f"{sample_id} ({date_str})"
        output_filename = f"{final_name}.docx"

        pdf_bytes = pdf_file.read()
        gene_names_list, percent_identity_str, ratio_str = extract_from_pdf(pdf_bytes)

        if gene_names_list and percent_identity_str and ratio_str:
            docx_bytes = generate_docx_report(
                gene_names_list=gene_names_list,
                percent_identity_str=percent_identity_str,
                ratio_str=ratio_str,
                template_bytes=docx_template_file.read(),
                sample_id_text=final_name
            )

            st.download_button(
                label="Download Report (.docx)",
                data=docx_bytes,
                file_name=output_filename,
                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document"
            )
        else:
            st.error("Could not extract required fields from PDF.")
    else:
        st.warning("Please upload all three required files.")
