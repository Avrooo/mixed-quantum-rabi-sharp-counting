set -e
cd /home/claude/pra/docx
python3 prep.py
pandoc body.tex -f latex -t docx -o body_raw.docx
pandoc abstract.tex -f latex -t docx -o abs_raw.docx
python3 assemble.py >/dev/null
cd pkg && rm -f ../out.docx && zip -XrDq ../out.docx . && cd ..
python3 /root/.claude/skills/synced/*/docx/scripts/office/validate.py out.docx 2>&1 | tail -2
python3 /root/.claude/skills/synced/*/docx/scripts/office/soffice.py --headless --convert-to pdf out.docx >/dev/null 2>&1
pdfinfo out.pdf | grep Pages
