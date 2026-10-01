import subprocess
import sys
from html import unescape
from pathlib import Path

from lxml import etree

from packages.contracts import parse_json

ROOT = Path(__file__).resolve().parents[2]
VERSION = "tpd-text-extractor-0.1.1"


def text_of(node):
    return " ".join(" ".join(node.itertext()).split())


def extract_xml(raw, max_segments):
    parser = etree.XMLParser(resolve_entities=False, no_network=True, load_dtd=False, huge_tree=False)
    root = etree.fromstring(raw, parser)
    if etree.QName(root).localname != "article":
        raise ValueError("XML_IS_NOT_JATS_ARTICLE")
    tree = root.getroottree()
    if tree.docinfo.internalDTD and tree.docinfo.internalDTD.entities():
        raise ValueError("XML_ENTITY_DECLARATIONS_UNSUPPORTED")
    nodes = root.xpath("//*[local-name()='article-title' or local-name()='p' or local-name()='table-wrap' or local-name()='fig']")
    chunks, warnings = [], []
    for node in nodes:
        name = etree.QName(node).localname
        ancestors = {etree.QName(a).localname for a in node.iterancestors() if isinstance(a.tag, str)}
        if "ref-list" in ancestors:
            continue
        if name == "p" and ("table-wrap" in ancestors or "fig" in ancestors or "ref-list" in ancestors):
            continue
        value = text_of(node)
        if not value:
            continue
        if len(chunks) >= max_segments:
            warnings.append("SEGMENT_LIMIT")
            break
        if len(value) > 12000:
            value = value[:12000]
            warnings.append("SEGMENT_TEXT_TRUNCATED")
        if name == "fig":
            warnings.append("FIGURES_NOT_INTERPRETED")
        if name == "table-wrap":
            warnings.append("TABLE_LAYOUT_NOT_INTERPRETED")
        chunks.append({"text": value, "xpath": tree.getpath(node), "content_kind": "figure_caption" if name == "fig" else "table" if name == "table-wrap" else "text",
                       "table_id": node.get("id") if name == "table-wrap" else None, "figure_id": node.get("id") if name == "fig" else None})
    return chunks, list(dict.fromkeys(warnings))


def extract_pdf(raw, max_segments):
    completed = subprocess.run([sys.executable, "-m", "packages.literature.pdf_worker", str(max_segments)],
        input=raw, capture_output=True, timeout=22, cwd=ROOT, check=True)
    result = parse_json(completed.stdout)
    warnings = ["PDF_LAYOUT_AND_IMAGES_NOT_INTERPRETED"]
    if result["truncated"]:
        warnings.append("PDF_PAGE_OR_TEXT_LIMIT")
    if not result["chunks"]:
        warnings.append("PDF_NO_EXTRACTABLE_TEXT")
    return result["chunks"], warnings


def strip_markup(text):
    root = etree.fromstring(("<div>" + unescape(text) + "</div>").encode("utf-8"),
        etree.HTMLParser(no_network=True, encoding="utf-8"))
    return text_of(root) if root is not None else text
