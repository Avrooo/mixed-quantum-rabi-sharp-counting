# Manuscript
`Sharp_counting_PRA_Word_template.docx` is the **sole authoritative manuscript** (Word, two-column APS layout, equations as native Word equation objects).

`Sharp_counting_PRA_equation_render.pdf` is **not** a Word export. It is a single-column rendering of the .docx content (pandoc -> MathML -> headless Chromium) produced only to check that every equation, figure, table and reference of the .docx is well-formed and legible. Layout, pagination and fonts differ from Word. The submission PDF must be exported from Word itself ("Save as PDF").

`build_tools_optional/` holds the scripts that originally generated the .docx from LaTeX sources. They are **not needed** to reproduce any result in this repository and the LaTeX sources are not shipped.
