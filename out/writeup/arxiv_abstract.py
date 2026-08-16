"""Emit the abstract as plain text for the arXiv web form.

arXiv's abstract field is plain text: LaTeX markup is not rendered, so macros
have to be resolved by hand rather than pasted.
"""
import pathlib
import re

src = pathlib.Path(__file__).with_name("manuscript.tex").read_text(encoding="utf8")
a = src.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]

subs = [
    (r"\\SIrange\{([^}]*)\}\{([^}]*)\}\{\\percent\}", r"\1-\2%"),
    (r"\\SI\{([^}]*)\}\{\\percent\}", r"\1%"),
    (r"\\SI\{([^}]*)\}\{[^}]*\}", r"\1"),
    (r"\\num\{([^}]*)\}", r"\1"),
    (r"\\nb", "nb"),
    (r"\\emph\{([^}]*)\}", r"\1"),
    (r"\\textbf\{([^}]*)\}", r"\1"),
    (r"\\Delta", "Delta "),
    (r"\\times", "x"),
    (r"\\leq", "<="),
    (r"\\cite\{[^}]*\}", ""),
    (r"\\[a-zA-Z]+", ""),
]
for pat, rep in subs:
    a = re.sub(pat, rep, a)

a = (a.replace("$", "").replace("{", "").replace("}", "")
      .replace("--", "-").replace("~", " "))
a = re.sub(r"[ \t]+", " ", a)
a = re.sub(r"\n\s*\n\s*", "\n\n", a).strip()
print(a)
print()
print(f"[{len(a)} characters; arXiv abstract limit is 1920]")
