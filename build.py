"""Build the multi-page handbook with Pandoc 3. Supports the archived and public layouts."""
from html import escape
from pathlib import Path
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "chapters" if (ROOT / "chapters").is_dir() else ROOT
OUT = ROOT / "docs" if SOURCE != ROOT else ROOT
TEMPLATE = ((ROOT / "web" if SOURCE != ROOT else ROOT) / "layout.html").read_text()
BOOK = "从 OLS 到因果识别"
chapters = []
index = []
paths = sorted(SOURCE.glob("[0-9][0-9]-*.md"))
for n, path in enumerate(paths, 1):
    title, content = path.read_text().split("\n", 1)
    title = title.removeprefix("# ")
    url = path.with_suffix(".html").name
    pieces = re.split(r"(?m)^## (.+)$", content)
    source = [pieces[0]]
    sections = []
    for k, (heading, body) in enumerate(zip(pieces[1::2], pieces[2::2]), 1):
        anchor = f"c{n}-s{k}"
        def section_link(match):
            chapter, section = match.groups()
            dest = paths[int(chapter)-1].with_suffix(".html").name
            return f"[{match.group()}]({dest}#c{chapter}-s{section})"
        linked_body = re.sub(r"第 (\d+)\.(\d+) 节", section_link, body)
        source.append(f"## {heading} {{#{anchor}}}\n{linked_body}")
        sections.append((heading, anchor))
        index.append([heading, url+"#"+anchor, re.sub(r"\s+", " ", body).lower()])
    result = subprocess.run(["pandoc", "-f", "markdown-smart", "-t", "html5", "--mathml", "--wrap=none"],
                            input="\n".join(source), text=True, capture_output=True, check=True)
    if result.stderr.strip():
        raise ValueError(result.stderr)
    chapters.append((title, url, sections, result.stdout))
OUT.mkdir(exist_ok=True)
def link(url, title):
    return f'<a href="{escape(url, quote=True)}">{escape(title)}</a>'

def render(file, title, active, body):
    nav = []
    for i, (name, url, *_) in enumerate(chapters, 1):
        current = ' aria-current="page"' if active == i else ''
        nav.append(f'<a class="chapterlink" href="{url}"{current}><b>{i:02}</b>{escape(name)}</a>')
    html = TEMPLATE.replace("@@TITLE@@", escape(title)).replace("@@NAV@@", "".join(nav)).replace("@@BODY@@", body)
    (OUT / file).write_text(html)

cards = []
for n, (name, url, sections, fragment) in enumerate(chapters, 1):
    toc = "".join("<li>"+link("#"+anchor, title)+"</li>" for title, anchor in sections)
    previous = link(chapters[n-2][1], "← "+chapters[n-2][0]) if n>1 else link("index.html", "← 全书目录")
    following = link(chapters[n][1], chapters[n][0]+" →") if n<len(chapters) else link("index.html", "全书目录 →")
    label = "第"+"一二三四五六七"[n-1]+"章"
    heading = f'<header id="c{n}"><a href="index.html">全书目录</a><span class="crumb"> / {label}</span><h1>{escape(name)}</h1></header>'
    render(url, name+" · "+BOOK, n, heading+f'<details class="contents"><summary>本章目录</summary><ul>{toc}</ul></details>'
           +'<article>'+fragment+'</article><nav class="pager" aria-label="章节翻页">'+previous+following+'</nav>')
    case = sections[-1][0].split(" ", 1)[1]
    cards.append(f'<section class="card"><span class="eyebrow">{label}</span><h2>{link(url,name)}</h2>'
                 +f'<p>{escape(case)}</p><details><summary>本章 {len(sections)} 节</summary><ul>'
                 +"".join("<li>"+link(url+"#"+anchor, title)+"</li>" for title, anchor in sections)+"</ul></details></section>")
intro = '<header><span class="eyebrow">计量经济学讲义</span><h1>'+BOOK+'</h1><p>概率与统计推断、线性投影和政策评估的识别条件。</p></header>'
render("index.html", BOOK+" · 全书目录", 0, intro+'<div class="cards">'+''.join(cards)+'</div>')
routes = [c[1] for c in chapters]
data = json.dumps(index, ensure_ascii=False)
(OUT / "search.js").write_text("const bookSections = "+data+";\nconst chapterFiles = "+json.dumps(routes)+";\n"+"""
const legacy = location.hash.match(/^#c([1-7])(?:-s[0-9]+)?$/);
if (legacy && !document.querySelector('article')) {
  location.replace(chapterFiles[Number(legacy[1])-1]+location.hash);
}
const searchBox = document.querySelector('#search');
searchBox.addEventListener('input', () => {
  const query = searchBox.value.trim().toLowerCase();
  const matches = query ? bookSections.filter(([title, url, body]) => (title.toLowerCase()+' '+body).includes(query)) : [];
  const list = document.querySelector('#results');
  list.replaceChildren();
  document.querySelector('#status').textContent = query ? '匹配 '+matches.length+' 节' : '';
  for (const [title, url] of matches) {
    const row = document.createElement('li');
    const a = document.createElement('a');
    a.href = url;
    a.textContent = title;
    row.append(a);
    list.append(row);
  }
});
""")
(OUT / ".nojekyll").touch()
