import streamlit as st
from io import BytesIO
from docx import Document
import PyPDF2
import re
from datetime import datetime

# ============================================================
# Helper function: extract data from PDF
# ============================================================
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

# ============================================================
# Helper function: generate DOCX
# ============================================================
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

    # Determine prognosis
    if any("3-21" in g for g in gene_names_list):
        prognosis_text = (
            "BAD\n"
            "Note: CLL expressing the IGHV 3-21 variable region gene segment "
            "have a poorer prognosis regardless of IGHV mutation status."
        )
        prognosis_single_line = "BAD"
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

    # Load DOCX template from bytes
    doc = Document(BytesIO(template_bytes))

    # Replace Sample ID in header
    if sample_id_text:
        for section in doc.sections:
            header = section.header
            for para in header.paragraphs:
                if "Sample ID" in para.text:
                    para.text = f"Sample ID: {sample_id_text}"

    # Replace Sample ID in body paragraphs
    if sample_id_text:
        for para in doc.paragraphs:
            if "Sample ID" in para.text:
                para.text = f"Sample ID: {sample_id_text}"

    # Replace Sample ID in tables
    if sample_id_text:
        for table in doc.tables:
            for row in table.rows:
                for cell in row.cells:
                    if "Sample ID" in cell.text:
                        cell.text = f"Sample ID: {sample_id_text}"

    # Replace mutation and prognosis in paragraphs
    for para in doc.paragraphs:
        if "Mutation detected:" in para.text:
            para.text = f"Mutation detected: {mutation_percent_str}"

        if "Clinical Prognosis" in para.text:
            idx = doc.paragraphs.index(para)
            if idx + 1 < len(doc.paragraphs):
                doc.paragraphs[idx + 1].text = prognosis_text

    # Replace mutation and prognosis in tables
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                if "Mutation detected:" in cell.text:
                    value_cell = row.cells[1]
                    value_cell.text = mutation_percent_str
                if "Clinical Prognosis:" in cell.text:
                    value_cell = row.cells[1]
                    value_cell.text = prognosis_single_line

    # Fill the Closest Allele table
    for table in doc.tables:
        first_row = [cell.text.strip() for cell in table.rows[0].cells]
        if (
            "Name" in first_row
            and "Percentage Identity Detected" in first_row
        ):
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

    # Save to BytesIO
    output = BytesIO()
    doc.save(output)
    output.seek(0)

    return output

# ============================================================
# STREAMLIT APP UI
# ============================================================
st.set_page_config(page_title="IGHV Report Generator", layout="centered")
st.title("IGHV Report Generator")

# ADD YOUR TAGLINE BELOW THE TITLE
st.caption("*A Wobble Base Bioresearch proprietary software*")

st.write("Please upload:")
pdf_file = st.file_uploader("PDF file (IgBLAST results)", type="pdf")
docx_template_file = st.file_uploader("Word Template (.docx)", type="docx")
ab1_file = st.file_uploader(".ab1 Sequence File", type="ab1")

if st.button("Generate Report"):
    if pdf_file and docx_template_file and ab1_file:
        # Extract sample ID
        ab1_filename = ab1_file.name
        sample_id = ab1_filename.split("(")[0].strip()
        date_str = datetime.now().strftime("%d-%m-%Y")
        final_name = f"{sample_id} ({date_str})"
        output_filename = f"{final_name}.docx"

        st.info(f"Sample ID: {sample_id}")
        st.info(f"Output file name: {output_filename}")

        # Process PDF
        pdf_bytes = pdf_file.read()
        gene_names_list, percent_identity_str, ratio_str = extract_from_pdf(pdf_bytes)

        st.write("Gene Names:", gene_names_list)
        st.write("Percent Identity:", percent_identity_str)
        st.write("Ratio:", ratio_str)

        if gene_names_list and percent_identity_str and ratio_str:
            # Generate DOCX
            docx_bytes = generate_docx_report(
                gene_names_list= gene_names_list,
                percent_identity_str=percent_identity_str,
                ratio_str=ratio_str,
                template_bytes=docx_template_file.read(),
                sample_id_text=final_name
            )

            # Provide download
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
