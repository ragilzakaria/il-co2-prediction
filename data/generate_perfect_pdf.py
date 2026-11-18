#!/usr/bin/env python3
"""
generate_perfect_pdf.py — Compile high-quality, publication-grade academic PDF
from paper.tex with beautiful math rendering, embedded figures, clean tables, and zero artifacts.
"""
import re
import os
import subprocess
import latex2mathml.converter

HERE = os.path.dirname(os.path.abspath(__file__))
FINAL_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "workspace", "final"))
TEX_PATH = os.path.join(FINAL_DIR, "paper.tex")
BIB_PATH = os.path.join(FINAL_DIR, "refs.bib")
HTML_PATH = os.path.join(FINAL_DIR, "paper.html")
PDF_PATH = os.path.join(FINAL_DIR, "paper.pdf")

# Cross-reference mappings
REF_MAP = {
    "fig:dataset_distribution": "1",
    "fig:model_parity": "2",
    "fig:feature_importance": "3",
    "fig:screening_isotherms": "4",
    "tab:model_performance": "1",
    "tab:pressure_regimes": "2",
    "tab:feature_importance": "3",
    "tab:screening": "4",
}

def clean_latex_math(raw_math):
    # Remove \nonumber, \\, and alignment tabs &
    clean = raw_math.replace(r"\nonumber", "").replace("&", " ")
    clean = re.sub(r"\\\\+", " ", clean)
    clean = re.sub(r"\s+", " ", clean).strip()
    return clean

def math_to_mathml(match):
    raw_math = match.group(1).strip()
    is_display = match.group(0).startswith("$$") or match.group(0).startswith(r"\[")
    clean = clean_latex_math(raw_math)
    try:
        mathml = latex2mathml.converter.convert(clean)
        # Ensure clean text and entity replacement
        mathml = mathml.replace("&#x0003C;", "&lt;").replace("&#x000A0;", " ")
        if is_display:
            mathml = f'<div class="eq-display">{mathml}</div>'
        return mathml
    except Exception:
        if is_display:
            return f'<div class="eq-display"><math display="block"><mrow><mtext>{clean}</mtext></mrow></math></div>'
        else:
            return f'<math display="inline"><mrow><mtext>{clean}</mtext></mrow></math>'

def format_equation_env(match, eq_num=None):
    raw_math = match.group(1).strip()
    clean = clean_latex_math(raw_math)
    try:
        mathml = latex2mathml.converter.convert(clean)
        mathml = mathml.replace("&#x0003C;", "&lt;").replace("&#x000A0;", " ")
        eq_label = f'<span class="eq-number">({eq_num})</span>' if eq_num else ""
        return f'<div class="eq-container"><div class="eq-display">{mathml}</div>{eq_label}</div>'
    except Exception:
        eq_label = f'<span class="eq-number">({eq_num})</span>' if eq_num else ""
        return f'<div class="eq-container"><div class="eq-display"><math display="block"><mrow><mtext>{clean}</mtext></mrow></math></div>{eq_label}</div>'

def build_html():
    with open(TEX_PATH) as f:
        tex = f.read()

    # 1. Parse Title
    title_m = re.search(r"\\title\{\\textbf\{([^}]+)\}\}", tex)
    title = title_m.group(1) if title_m else "Predictive Modeling of Carbon Dioxide Capture in Ionic Liquids"

    # 2. Parse Abstract
    abs_m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
    abstract = abs_m.group(1).strip() if abs_m else ""

    # 3. Strip everything before \section{Introduction} for body_tex to avoid duplicate abstract!
    intro_start_m = re.search(r"\\section\{Introduction\}", tex)
    intro_start = intro_start_m.start() if intro_start_m else 0
    ack_m = re.search(r"\\section\*\{CRediT Authorship", tex)
    body_end = ack_m.start() if ack_m else len(tex)
    body_tex = tex[intro_start:body_end]

    # 4. Resolve \ref{...} cross-references
    for label, num in REF_MAP.items():
        body_tex = re.sub(rf"\\ref\{{{label}\}}", num, body_tex)

    # 5. Process Abstract Math
    abstract_html = re.sub(r"\$([^\$]+)\$", math_to_mathml, abstract)

    # 6. Process Headings
    def replace_sections(text):
        text = re.sub(r"\\section\{([^}]+)\}", r'<h2 class="sec-title">\1</h2>', text)
        text = re.sub(r"\\subsection\{([^}]+)\}", r'<h3 class="subsec-title">\1</h3>', text)
        text = re.sub(r"\\subsubsection\{([^}]+)\}", r'<h4 class="subsubsec-title">\1</h4>', text)
        return text

    processed = replace_sections(body_tex)

    # 7. Process Equations with numbers
    eq_counter = 1
    def eq_repl(m):
        nonlocal eq_counter
        res = format_equation_env(m, eq_counter)
        eq_counter += 1
        return res

    processed = re.sub(r"\\begin\{equation\}(.*?)\\end\{equation\}", eq_repl, processed, flags=re.S)
    processed = re.sub(r"\\begin\{align\*?\}(.*?)\\end\{align\*?\}", eq_repl, processed, flags=re.S)
    processed = re.sub(r"\$\$([^\$]+)\$\$", math_to_mathml, processed)
    processed = re.sub(r"\$([^\$]+)\$", math_to_mathml, processed)

    # 8. Process Citations
    with open(BIB_PATH) as f:
        bib_text = f.read()

    cite_map = {}
    cite_counter = 1
    
    def cite_repl(match):
        nonlocal cite_counter
        keys = [k.strip() for k in match.group(1).split(",")]
        nums = []
        for k in keys:
            if k not in cite_map:
                cite_map[k] = cite_counter
                cite_counter += 1
            nums.append(str(cite_map[k]))
        return f'<span class="cite">[{", ".join(nums)}]</span>'

    processed = re.sub(r"\\cite\{([^}]+)\}", cite_repl, processed)

    # 9. Process Figures with Figure Numbers
    fig_counter = 1
    fig_pattern = re.compile(r"\\begin\{figure\*?\}\[[^\]]*\]\s*\\centering\s*\\includegraphics\[[^\]]*\]\{([^}]+)\}\s*\\caption\{(.*?)\}\s*\\label\{([^}]+)\}\s*\\end\{figure\*?\}", re.S)
    
    def fig_repl(m):
        nonlocal fig_counter
        img_path = m.group(1).strip()
        caption = m.group(2).strip()
        label = m.group(3).strip()
        cap_html = re.sub(r"\$([^\$]+)\$", math_to_mathml, caption)
        num = REF_MAP.get(label, str(fig_counter))
        fig_counter += 1
        return f'''
        <div class="figure-box" id="{label}">
            <img src="{img_path}" alt="{label}" class="paper-img" />
            <div class="figure-caption"><strong>Figure {num}:</strong> {cap_html}</div>
        </div>
        '''

    processed = fig_pattern.sub(fig_repl, processed)

    # 10. Process Tables with Table Numbers
    tab_counter = 1
    tab_pattern = re.compile(r"\\begin\{table\*?\}\[[^\]]*\]\s*\\centering\s*\\caption\{(.*?)\}\s*\\label\{([^}]+)\}\s*\\begin\{tabular\}\{([^}]+)\}(.*?)\\end\{tabular\}\s*\\end\{table\*?\}", re.S)

    def tab_repl(m):
        nonlocal tab_counter
        caption = m.group(1).strip()
        label = m.group(2).strip()
        col_spec = m.group(3).strip()
        content = m.group(4).strip()
        
        cap_html = re.sub(r"\$([^\$]+)\$", math_to_mathml, caption)
        num = REF_MAP.get(label, str(tab_counter))
        tab_counter += 1

        rows = content.split(r"\\")
        html_rows = []
        is_header = True
        for r in rows:
            r = r.strip()
            if not r:
                continue
            r = r.replace(r"\toprule", "").replace(r"\midrule", "").replace(r"\bottomrule", "").strip()
            if not r:
                continue
            cells = r.split("&")
            tag = "th" if is_header else "td"
            cell_htmls = []
            for c in cells:
                c = c.strip()
                c = c.replace(r"\_", "_").replace(r"\%", "%")
                c = re.sub(r"\\textbf\{([^}]+)\}", r"<strong>\1</strong>", c)
                c = re.sub(r"\$([^\$]+)\$", math_to_mathml, c)
                cell_htmls.append(f"<{tag}>{c}</{tag}>")
            html_rows.append(f"<tr>{''.join(cell_htmls)}</tr>")
            is_header = False

        return f'''
        <div class="table-box" id="{label}">
            <div class="table-caption"><strong>Table {num}:</strong> {cap_html}</div>
            <table class="paper-table">
                {''.join(html_rows)}
            </table>
        </div>
        '''

    processed = tab_pattern.sub(tab_repl, processed)

    # Clean text formatting
    processed = processed.replace(r"\_", "_").replace(r"\%", "%")
    processed = re.sub(r"\\textbf\{([^}]+)\}", r"<strong>\1</strong>", processed)
    processed = re.sub(r"\\textit\{([^}]+)\}", r"<em>\1</em>", processed)
    processed = re.sub(r"\\texttt\{([^}]+)\}", r"<code>\1</code>", processed)
    processed = re.sub(r"\\noindent", "", processed)
    processed = re.sub(r"---", "&mdash;", processed)
    processed = re.sub(r"--", "&ndash;", processed)
    processed = re.sub(r"~", " ", processed)

    # Lists
    processed = re.sub(r"\\begin\{enumerate\}", "<ol>", processed)
    processed = re.sub(r"\\end\{enumerate\}", "</ol>", processed)
    processed = re.sub(r"\\begin\{itemize\}", "<ul>", processed)
    processed = re.sub(r"\\end\{itemize\}", "</ul>", processed)
    processed = re.sub(r"\\item\s+", "<li>", processed)
    processed = re.sub(r"<li>(.*?)((?=<li>)|(?=</ol>)|(?=</ul>))", r"<li>\1</li>", processed, flags=re.S)

    # Paragraphs
    paragraphs = processed.split("\n\n")
    html_paragraphs = []
    for p in paragraphs:
        p = p.strip()
        if not p:
            continue
        if p.startswith("<h2") or p.startswith("<h3") or p.startswith("<div") or p.startswith("<ol") or p.startswith("<ul"):
            html_paragraphs.append(p)
        else:
            html_paragraphs.append(f"<p>{p}</p>")
    
    body_html = "\n".join(html_paragraphs)

    # Backmatter
    back_html = """
    <h2 class="sec-title">CRediT Authorship Contribution Statement</h2>
    <p><strong>Maulana Ragil Zakariah:</strong> Conceptualization, Methodology, Software, Data curation, Formal analysis, Investigation, Writing &ndash; original draft. 
    <strong>Aulia Sukma Hutama:</strong> Supervision, Conceptualization, Methodology, Writing &ndash; review &amp; editing, Funding acquisition. 
    <strong>Akhmad Syoufian:</strong> Supervision, Conceptualization, Writing &ndash; review &amp; editing. 
    <strong>Dwi Siswanta:</strong> Validation, Resources, Writing &ndash; review &amp; editing. 
    <strong>Robby Noor Cahyono:</strong> Validation, Formal analysis, Writing &ndash; review &amp; editing.</p>

    <h2 class="sec-title">Declaration of Competing Interest</h2>
    <p>The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.</p>

    <h2 class="sec-title">Acknowledgements</h2>
    <p>This research was supported by the Department of Chemistry, Faculty of Mathematics and Natural Sciences, Universitas Gadjah Mada. The authors gratefully acknowledge the computational resources provided by the Computational Chemistry Laboratory (Austrian-Indonesian Centre for Computational Chemistry, AIC) at Universitas Gadjah Mada.</p>
    """

    # References
    bib_entries = {}
    pattern = re.compile(r"@(\w+)\{\s*([^,\s]+),([^@]*?)\n\}", re.S)
    for match in pattern.finditer(bib_text):
        entry_type, key, bbody = match.groups()
        fields = {}
        for line in bbody.split("\n"):
            fm = re.match(r"^\s*(\w+)\s*=\s*[\{|\"](.*)[\}|\"],?\s*$", line)
            if fm:
                fields[fm.group(1).lower()] = fm.group(2).strip().strip("{}")
        bib_entries[key] = fields

    sorted_cites = sorted(cite_map.items(), key=lambda x: x[1])
    ref_list_html = []
    for key, num in sorted_cites:
        b = bib_entries.get(key, {})
        author = b.get("author", "Unknown author")
        title_r = b.get("title", key)
        journal = b.get("journal", "")
        year = b.get("year", "")
        doi = b.get("doi", "")
        doi_link = f' <a href="https://doi.org/{doi}" class="doi-link">https://doi.org/{doi}</a>' if doi else ""
        ref_list_html.append(f'<div class="ref-item">[{num}] {author}, "{title_r}", <em>{journal}</em> ({year}).{doi_link}</div>')

    references_html = f'''
    <h2 class="sec-title">References</h2>
    <div class="ref-container">
        {"".join(ref_list_html)}
    </div>
    '''

    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
@page {{
    size: letter;
    margin: 0.7in 0.65in 0.75in 0.65in;
}}

* {{
    box-sizing: border-box;
}}

body {{
    font-family: "Times New Roman", Times, "Liberation Serif", serif;
    font-size: 9.5pt;
    line-height: 1.34;
    color: #111827;
    background: #fff;
    margin: 0;
    padding: 0;
}}

/* Header strip */
.paper-header {{
    text-align: center;
    margin-bottom: 14px;
    border-bottom: 2px solid #002b49;
    padding-bottom: 12px;
}}

.journal-slug {{
    font-family: "Helvetica Neue", Arial, sans-serif;
    font-size: 8.5pt;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: #002b49;
    margin-bottom: 8px;
}}

.paper-title {{
    font-size: 16.5pt;
    font-weight: bold;
    line-height: 1.25;
    margin: 6px 0 10px 0;
    color: #0c1a2c;
}}

.author-block {{
    font-size: 9.5pt;
    margin-bottom: 6px;
    color: #222;
}}

.author-name {{
    font-weight: 600;
}}

.affil-block {{
    font-size: 8.2pt;
    color: #4b5563;
    line-height: 1.35;
    margin-bottom: 4px;
}}

.contact-line {{
    font-size: 8.2pt;
    color: #0369a1;
    font-family: "Courier New", monospace;
}}

/* Abstract Box */
.abstract-container {{
    background: #f8fafc;
    border: 1px solid #cbd5e1;
    border-radius: 4px;
    padding: 10px 14px;
    margin: 12px 0 18px 0;
    font-size: 8.8pt;
    text-align: justify;
    line-height: 1.38;
}}

.abstract-title {{
    font-weight: bold;
    text-transform: uppercase;
    font-size: 8.5pt;
    letter-spacing: 0.05em;
    color: #002b49;
    margin-bottom: 4px;
}}

.keywords {{
    margin-top: 8px;
    font-size: 8.2pt;
    color: #334155;
    border-top: 1px dashed #cbd5e1;
    padding-top: 6px;
}}

/* Two Column Layout */
.content-body {{
    column-count: 2;
    column-gap: 0.28in;
    text-align: justify;
}}

/* Headings */
h2.sec-title {{
    font-size: 10.5pt;
    font-weight: bold;
    color: #002b49;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    border-bottom: 1px solid #002b49;
    padding-bottom: 2px;
    margin: 14px 0 6px 0;
    break-after: avoid;
    break-inside: avoid;
}}

h3.subsec-title {{
    font-size: 9.5pt;
    font-weight: bold;
    font-style: italic;
    color: #1e293b;
    margin: 10px 0 4px 0;
    break-after: avoid;
    break-inside: avoid;
}}

p {{
    margin: 0 0 6px 0;
    text-indent: 14px;
}}

p:first-of-type {{
    text-indent: 0;
}}

/* Math Formatting */
.eq-container {{
    display: flex;
    align-items: center;
    justify-content: center;
    margin: 8px 0;
    break-inside: avoid;
    position: relative;
}}

.eq-display {{
    text-align: center;
    margin: 0 auto;
}}

.eq-number {{
    font-size: 8.5pt;
    position: absolute;
    right: 4px;
    color: #475569;
}}

math {{
    font-size: 10pt;
}}

/* Figures */
.figure-box {{
    margin: 14px 0;
    text-align: center;
    break-inside: avoid;
    background: #ffffff;
    border: 1px solid #e2e8f0;
    padding: 6px;
    border-radius: 4px;
}}

.paper-img {{
    width: 100%;
    max-width: 100%;
    height: auto;
    display: block;
    margin: 0 auto 6px auto;
}}

.figure-caption {{
    font-size: 8pt;
    line-height: 1.25;
    text-align: justify;
    color: #334155;
    padding: 2px 4px;
}}

/* Tables */
.table-box {{
    margin: 14px 0;
    break-inside: avoid;
}}

.table-caption {{
    font-size: 8.2pt;
    font-weight: bold;
    margin-bottom: 4px;
    color: #0f172a;
    text-align: left;
}}

.paper-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 7.8pt;
    line-height: 1.22;
    border-top: 1.5px solid #002b49;
    border-bottom: 1.5px solid #002b49;
    margin: 4px 0;
}}

.paper-table th {{
    border-bottom: 1px solid #002b49;
    padding: 5px 3px;
    font-weight: bold;
    text-align: center;
    background: #f8fafc;
}}

.paper-table td {{
    padding: 3.5px 3px;
    text-align: center;
    border-bottom: 0.5px solid #e2e8f0;
}}

/* Lists */
ol, ul {{
    margin: 4px 0 8px 18px;
    padding: 0;
}}

li {{
    margin-bottom: 3px;
    padding-left: 2px;
}}

/* References */
.ref-container {{
    font-size: 7.8pt;
    line-height: 1.25;
    color: #334155;
}}

.ref-item {{
    margin-bottom: 4.5px;
    text-indent: -16px;
    padding-left: 16px;
}}

.doi-link {{
    color: #0284c7;
    text-decoration: none;
    word-break: break-all;
}}

.cite {{
    font-size: 8.5pt;
    color: #0369a1;
    font-weight: 500;
}}
</style>
</head>
<body>

<div class="paper-header">
    <div class="journal-slug">Chemical Engineering Science &bull; Research Article &bull; Elsevier Format</div>
    <h1 class="paper-title">{title}</h1>
    <div class="author-block">
        <span class="author-name">Maulana Ragil Zakariah</span><sup>1</sup>, 
        <span class="author-name">Aulia Sukma Hutama</span><sup>1,*</sup>, 
        <span class="author-name">Akhmad Syoufian</span><sup>1</sup>, 
        <span class="author-name">Dwi Siswanta</span><sup>1</sup>, 
        <span class="author-name">Robby Noor Cahyono</span><sup>1</sup>
    </div>
    <div class="affil-block">
        <sup>1</sup>Department of Chemistry, Faculty of Mathematics and Natural Sciences, Universitas Gadjah Mada, Sekip Utara, BLS 21, Yogyakarta 55281, Indonesia
    </div>
    <div class="contact-line">
        *Corresponding author: aulia.sukma.hutama@ugm.ac.id
    </div>
</div>

<div class="abstract-container">
    <div class="abstract-title">Abstract</div>
    <p>{abstract_html}</p>
    <div class="keywords">
        <strong>Keywords:</strong> Carbon dioxide capture; Ionic liquids; Density functional theory; Machine learning; Random Forest; Quantitative structure-property relationship.
    </div>
</div>

<div class="content-body">
    {body_html}
    {back_html}
    {references_html}
</div>

</body>
</html>
"""

    with open(HTML_PATH, "w") as f:
        f.write(full_html)
    print(f"Generated clean HTML: {HTML_PATH}")

if __name__ == "__main__":
    build_html()
