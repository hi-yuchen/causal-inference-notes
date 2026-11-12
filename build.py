"""Build a self-contained GitHub Pages site. Requires Python 3 and Pandoc 3."""
from html import escape
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parent
template = (root / "layout.html").read_text(encoding="utf-8")
sections, links = [], []
for n, path in enumerate(sorted(root.glob("[0-9][0-9]-*.md")), 1):
    lines = path.read_text(encoding="utf-8").splitlines()
    title = lines[0].removeprefix("# ")
    toc, paragraphs = [], []
    for line in lines[1:]:
        if line.startswith("## "):
            heading = line[3:]
            slug = f"c{n}-s{len(toc)+1}"
            toc.append(f'<a href="#{slug}">{escape(heading)}</a>')
            line = f'### {heading} {{#{slug}}}'
        paragraphs.append(line)
    fragment = subprocess.run(
        ["pandoc", "-f", "markdown-smart", "-t", "html5", "--mathml", "--wrap=none"],
        input="\n".join(paragraphs), capture_output=True, text=True, check=True)
    if fragment.stderr.strip():
        raise RuntimeError(fragment.stderr)
    tag = f"第{'一二三四五六'[n-1]}篇" if n < 7 else "终章"
    links.append(f'<a class="chapterlink" href="#c{n}"><b>{n:02}</b>{escape(title)}</a>')
    sections.append(f'<article class="chapter" id="c{n}"><span class="eyebrow">{tag}</span>'
                    f'<h2>{escape(title)}</h2><details class="contents"><summary>本篇目录</summary>'
                    + ''.join(toc) + '</details>' + fragment.stdout
                    + '<a class="back" href="#home">返回总目录 ↑</a></article>')
output = root
output.mkdir(exist_ok=True)
(output / "index.html").write_text(template.replace("@@NAV@@", "".join(links))
    .replace("@@BODY@@", "\n".join(sections)), encoding="utf-8")
(output / ".nojekyll").touch()
