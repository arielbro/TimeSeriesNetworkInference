#!/usr/bin/env python3
"""Bulk-export every published Boolean model from Cell Collective as truth tables.

Endpoints mirror the official ccapi client (Cell Collective API v2):
    list   : GET {BASE}/api/model/cards/research
             ?modelTypes=boolean&orderBy=recent&category=published&cards=N
    export : GET {BASE}/_api/model/export/{id}?version={v}&type=TT
             (note the `_api` prefix on export vs `api` on the card list)

Export-type tokens are UPPERCASE: SBML | TT | EXPR | MATRIX | GML.
`type=TT` returns a zip whose contents match the truth-table folder layout
(per-component CSVs, SPECIES_KEY.csv, external_components_ALL.txt).

Only the `requests` library is required.
"""
import os
import time
import zipfile
import requests

BASE      = os.environ.get("URL", "https://cellcollective.org")
EXPORT    = "TT"                              # swap for SBML / EXPR / MATRIX / GML
OUTDIR    = "cellcollective_models"
UNZIP     = True                              # also expand each zip into OUTDIR/{id}/
MAX_CARDS = 10000                             # upper bound on models to list in one call

session = requests.Session()
session.headers["User-Agent"] = "cc-bulk-export/1.0"


def list_boolean_models():
    """Return [(model_id, version), ...] for every published Boolean model."""
    r = session.get(
        f"{BASE}/api/model/cards/research",
        params=[
            ("modelTypes", "boolean"),
            ("orderBy",    "recent"),
            ("category",   "published"),
            ("cards",      MAX_CARDS),
        ],
        timeout=60,
    )
    r.raise_for_status()
    cards = r.json()

    # Tolerate either a bare list or a wrapper object.
    if isinstance(cards, dict):
        cards = (cards.get("cards") or cards.get("content")
                 or next((v for v in cards.values() if isinstance(v, list)), []))

    out = []
    for c in cards:
        m = c.get("model", c)                       # card may nest fields under "model"
        vmap = m.get("modelVersionMap")
        if vmap:
            versions = list(vmap.keys())
        else:
            versions = [v["version"] for v in c.get("versions", [{"version": 1}])]
        out.append((m["id"], max(int(v) for v in versions)))   # latest version
    return out


def export(model_id, version):
    """Download one model's export to cc_export/{id}.zip (+ optional unzip)."""
    r = session.get(
        f"{BASE}/_api/model/export/{model_id}",
        params={"version": version, "type": EXPORT},
        stream=True, timeout=120,
    )
    r.raise_for_status()
    path = os.path.join(OUTDIR, f"{model_id}.zip")
    with open(path, "wb") as fh:
        for chunk in r.iter_content(chunk_size=8192):
            fh.write(chunk)
    if UNZIP:
        with zipfile.ZipFile(path) as z:
            z.extractall(os.path.join(OUTDIR, str(model_id)))
    return path


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    models = list_boolean_models()
    print(f"{len(models)} published Boolean models")
    for i, (mid, ver) in enumerate(models, 1):
        try:
            export(mid, ver)
            print(f"[{i}/{len(models)}] {mid} v{ver} ok")
        except Exception as e:
            print(f"[{i}/{len(models)}] {mid} v{ver} FAILED: {e}")
        time.sleep(0.3)                              # be gentle on the server


if __name__ == "__main__":
    main()