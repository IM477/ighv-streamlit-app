import streamlit as st
from io import BytesIO
from docx import Document
import PyPDF2
import re
from datetime import datetime
from Bio.Seq import Seq

# ============================================================
# SECTION 1: UGENE-style Consensus Generator from TXT file
# ============================================================

st.set_page_config(page_title="IGHV Tools", layout="centered")
st.title("Consensus Generator from Forward and Reverse Reads")

txt_file = st.file_uploader("Upload Forward/Reverse Read TXT file", type="txt", key="txt_input")
threshold = st.number_input("Set mismatch threshold (T)", min_value=0, value=10, step=1)

def parse_fasta_txt(file_content):
    content = file_content.decode("utf-8")
    entries = content.strip().split('>')
    seqs = {}
    for entry in entries:
        if not entry:
            continue
        lines = entry.strip().splitlines()
        header = lines[0].strip()
        sequence = ''.join(lines[1:]).replace(" ", "").upper()
        seqs[header] = sequence
    return seqs

def ugene_style_consensus(forward, reverse, threshold=10):
    s1 = forward.strip().replace("\n", "").upper()
    s2 = reverse.strip().replace("\n", "").upper()
    s2_rev = s2[::-1]
    s3 = str(Seq(s2_rev).complement())

    first_match = None
    last_match = None
    for i in range(min(len(s1), len(s3))):
        if s1[i] == s3[i]:
            if first_match is None:
                first_match = i
            last_match = i

    mismatch_count = 0
    a2 = ""
    if first_match is not None and last_match is not None:
        for i in range(first_match, last_match + 1):
            if s1[i] != s3[i]:
                mismatch_count += 1
        a2 = s1[first_match:last_match + 1]
    else:
        return None, s1, s2, s2_rev, s3, None, None, 0, None, None

    if mismatch_count > threshold:
        return None, s1, s2, s2_rev, s3, first_match, last_match, mismatch_count, None, None

    prefix = s1[:last_match + 1]
    suffix = s3[last_match + 1:]
    a3 = suffix
    consensus = prefix + suffix
    return consensus, s1, s2, s2_rev, s3, first_match, last_match, mismatch_count, a2, a3

if txt_file:
    seq_data = parse_fasta_txt(txt_file.read())

    forward = next((v for k, v in seq_data.items() if k.endswith("_F")), None)
    reverse = next((v for k, v in seq_data.items() if k.endswith("_R")), None)

    if forward and reverse:
        consensus, s1, s2, s2_rev, s3, first_match, last_match, mismatch_count, a2, a3 = ugene_style_consensus(forward, reverse, threshold)

        st.subheader("Inputs and Intermediates")
        st.text_area("Forward Read (s1)", s1, height=100)
        st.text_area("Reverse Read (s2)", s2, height=100)
        st.text_area("Reverse of s2 (s2_rev)", s2_rev, height=100)
        st.text_area("Reverse Complement (s3)", s3, height=100)

        if first_match is not None and last_match is not None:
            st.write(f"🔍 First Match Position: {first_match}")
            st.write(f"🔍 Last Match Position: {last_match}")
            st.write(f"❌ Mismatches: {mismatch_count}")
        else:
            st.error("❌ No matching region found between s1 and s3.")

        if consensus:
            st.subheader("✅ Consensus Output")
            st.text_area("Matching region (a2)", a2, height=100)
            st.text_area("Unmatched suffix from s3 (a3)", a3, height=100)
            st.text_area("Final Consensus Sequence", consensus, height=150)

            consensus_bytes = BytesIO()
            consensus_bytes.write(consensus.encode('utf-8'))
            consensus_bytes.seek(0)

            st.download_button(
                label="Download Consensus",
                data=consensus_bytes,
                file_name="consensus_output.txt",
                mime="text/plain"
            )
        elif first_match is not None:
            st.warning(f"⚠️ Mismatch count = {mismatch_count} exceeds threshold = {threshold}")
    else:
        st.error("❌ Could not find both IGHV_F and IGHV_R sequences in TXT.")

# ============================================================
# SECTION 2: IGHV Report Generator (original logic preserved)
# ============================================================

st.title("IGHV Report Generator")
st.caption("*A Wobble Base Bioresearch proprietary software*")

st.write("Please upload:")
pdf_file = st.file_uploader("PDF file (IgBLAST results)", type="pdf")
docx_template_file = st.file_uploader("Word Template (.docx)", type="docx")
ab1_file = st.file_uploader(".ab1 Sequence File", type="ab1")

def extract_from_pdf(pdf_bytes):
    pdf_reader = PyPDF2.PdfReader(BytesIO(pdf_bytes))
    text = ""
    for page in pdf_reader.pages:
        page_text = page.extract_text()
        text += "\n" + page_text

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
            raw_line = line.strip()
            if not raw_line:
                continue
            words = raw_line.split()
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

def generate_docx_report(gene_names_list, percent_identity_str, ratio_str, template_bytes, sample_id_text):
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
            prognosis_single_line = "BAD"
        elif 2.1 <= mutation_percent <= 3.0:
            prognosis_text = "Borderline with intermediate clinical course"
            prognosis_single_line = "Borderline with intermediate clinical course"
        else:
            prognosis_text = "GOOD"
            prognosis_single_line = "GOOD"

    doc = Document(BytesIO(template_bytes))

    if sample_id_text:
        for section in doc.sections:
            header = section.header
            for para in header.paragraphs:
                if "Sample ID" in para.text:
                    para.text = f"Sample ID: {sample_id_text}"
        for para in doc.paragraphs:
            if "Sample ID" in para.text:
                para.text = f"Sample ID: {sample_id_text}"
        for table in doc.tables:
            for row in table.rows:
                for cell in row.cells:
                    if "Sample ID" in cell.text:
                        cell.text = f"Sample ID: {sample_id_text}"

    for para in doc.paragraphs:
        if "Mutation detected:" in para.text:
            para.text = f"Mutation detected: {mutation_percent_str}"
        if "Clinical Prognosis" in para.text:
            idx = doc.paragraphs.index(para)
            if idx + 1 < len(doc.paragraphs):
                doc.paragraphs[idx + 1].text = prognosis_text

    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                if "Mutation detected:" in cell.text:
                    value_cell = row.cells[1]
                    value_cell.text = mutation_percent_str
                if "Clinical Prognosis:" in cell.text:
                    value_cell = row.cells[1]
                    value_cell.text = prognosis_single_line

    for table in doc.tables:
        first_row = [cell.text.strip() for cell in table.rows[0].cells]
        if "Name" in first_row and "Percentage Identity Detected" in first_row:
            if len(table.rows) > 1:
                for i in range(len(table.rows) - 1, 0, -1):
                    tbl = table._tbl
                    tr = table.rows[i]._tr
                    tbl.remove(tr)
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

if st.button("Generate Report"):
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
