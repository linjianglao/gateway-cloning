#!/usr/bin/env python3
"""Assemble the Gateway Bench sources into two single-file builds.

    index.html            standalone page for GitHub Pages / offline use
    dist/artifact.html    body-only build for claude.ai Artifacts
                          (the viewer supplies <!doctype>, <head> and a reset)

Usage:  python3 build.py
"""
import pathlib

ROOT = pathlib.Path(__file__).parent
SRC = ROOT / "src"
ORDER = ["core.js", "vectors.js", "lib.js", "lib2.js", "app.js"]

STANDALONE_HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<style>
:root{color-scheme:light;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}
body{margin:0;font-family:system-ui,-apple-system,Segoe UI,sans-serif;font-size:14px;background:#EEF1F4}
img{max-width:100%}
[hidden]{display:none!important}
</style>
"""


def script() -> str:
    parts = ["\n<script>\n"]
    for name in ORDER:
        parts.append("\n/* ==== %s ==== */\n" % name)
        parts.append((SRC / name).read_text())
    parts.append("\n</script>\n")
    return "".join(parts)


def main() -> None:
    head = (SRC / "main-head.html").read_text()
    body = head + script()

    # claude.ai artifact build: no <html>/<head>/<body> of its own
    dist = ROOT / "dist"
    dist.mkdir(exist_ok=True)
    (dist / "artifact.html").write_text(body)

    # standalone build: lift the page's own <title>/<link>/<style> into <head>
    cut = body.index("</style>") + len("</style>")
    page = (
        STANDALONE_HEAD + body[:cut] + "\n</head>\n<body>\n"
        + body[cut:] + "\n</body>\n</html>\n"
    )
    (ROOT / "index.html").write_text(page)

    print("index.html         %7d bytes" % len(page))
    print("dist/artifact.html %7d bytes" % len(body))


if __name__ == "__main__":
    main()
