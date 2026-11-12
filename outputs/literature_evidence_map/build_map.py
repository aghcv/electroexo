"""Build an auditable radial map from the supplied literature tracker.

Extraction uses bundled openpyxl read-only. Rendering uses Matplotlib and the
repository's manuscript styling. This never edits the source workbook/draft.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "figures" / "literature_evidence_map"
SOURCE = Path("/Users/aghorban/Library/CloudStorage/OneDrive-OldDominionUniversity/Papers/preparation/electro-exocytosis-comp-framework/nsPEF_EV_literature_evidence.xlsx")
DRAFT = Path("/Users/aghorban/code/manuscripts/electroexocytosis-computational-plos-computational-biology")
DATA = Path(__file__).with_name("extracted_records.json")
MODULE_SHORT = {
    1: "Pulse delivery and dosimetry",
    2: "Membrane and organelle electrodynamics",
    3: r"Ion transport, Ca$^{2+}$, ROS and ATP",
    4: r"Ca$^{2+}$-dependent remodeling and repair",
    5: "EV biogenesis and release",
    6: "Cargo, composition and potency",
    7: "Injury, debris and quality gate",
    8: "Manufacturing, isolation and QC",
}


def extract() -> None:
    import openpyxl

    workbook = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    sheet = workbook["Literature Tracker"]
    values = list(sheet.values)
    header = values[0]
    records = []
    within_table = collections.Counter()
    for excel_row, cells in enumerate(values[1:], start=2):
        if not any(x is not None for x in cells):
            continue
        raw = dict(zip(header, cells))
        table = raw["Table number"]
        within_table[table] += 1
        role = "Computational framework" if raw["Table type"] == "Computational framework" else "Biology/process evidence"
        slots = []
        for n in range(1, 6):
            citation = raw[f"Paper {n}"]
            if citation is None or not str(citation).strip():
                state = "Empty"
            elif len(str(citation).strip()) < 12:
                state = "Invalid placeholder"
            else:
                state = "Listed citation"
            slots.append({"number": n, "state": state, "raw": citation,
                          "cell": f"{openpyxl.utils.get_column_letter(11+n)}{excel_row}"})
        record = {
            "figure_id": f"{len(records)+1:02d}",
            "record_id": f"{table}.{within_table[table]:02d}",
            "source_excel_row": excel_row,
            "source_table": table,
            "source_table_type": raw["Table type"],
            "layer": int(raw["Layer"]),
            "module": raw["Module"],
            "submodule": raw["Submodule"],
            "role": role,
            "confidence": raw["Confidence"] if raw["Confidence"] else "Not assigned to model rows",
            "status": raw["Status"],
            "citation_count": sum(s["state"] == "Listed citation" for s in slots),
            "slots": slots,
            "model_form_or_phenomenon": raw["Biological phenomenon / Model form"],
            "evidence_focus": raw["Evidence to curate / Model focus"],
            "representative_relation": raw["Representative relation"],
            "inputs": raw["Inputs"],
            "outputs_or_relevance": raw["Outputs / Relevance"],
        }
        records.append(record)
    workbook.close()
    summary = {
        "records": len(records),
        "modules": len(set(r["layer"] for r in records)),
        "roles": dict(collections.Counter(r["role"] for r in records)),
        "biology_confidence": dict(collections.Counter(r["confidence"] for r in records if r["role"] != "Computational framework")),
        "statuses": dict(collections.Counter(r["status"] for r in records)),
        "occupied_citation_slots": sum(s["state"] != "Empty" for r in records for s in r["slots"]),
        "usable_citation_entries": sum(r["citation_count"] for r in records),
        "invalid_placeholders": [{"figure_id": r["figure_id"], "record_id": r["record_id"], "cell": s["cell"], "raw": s["raw"]} for r in records for s in r["slots"] if s["state"] == "Invalid placeholder"],
        "records_without_citations": [r["record_id"] for r in records if not r["citation_count"]],
        "module_summary": [{"layer": layer, "module": next(r["module"] for r in records if r["layer"] == layer),
                            "records": sum(r["layer"] == layer for r in records),
                            "usable_citation_entries": sum(r["citation_count"] for r in records if r["layer"] == layer)} for layer in range(1,9)],
    }
    # Expected-source checks detect a changed workbook rather than silently
    # publishing stale totals. These are extraction controls, not study scores.
    assert summary["records"] == 95
    assert summary["usable_citation_entries"] == 423
    assert summary["roles"] == {"Biology/process evidence": 47, "Computational framework": 48}
    assert summary["biology_confidence"] == {"Established": 23, "Supported": 19, "Hypothesized": 5}
    assert summary["records_without_citations"] == ["A3.01", "A3.02", "A3.03", "A3.04"]
    assert summary["invalid_placeholders"] == [{"figure_id": "51", "record_id": "A9.06", "cell": "N52", "raw": "s"}]
    payload = {"source": str(SOURCE), "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
               "manuscript": str(DRAFT / "main.tex"), "supplement": str(DRAFT / "supporting_information.tex"),
               "summary": summary, "records": records}
    DATA.parent.mkdir(parents=True, exist_ok=True)
    DATA.write_text(json.dumps(payload, ensure_ascii=False, indent=2))
    OUT.mkdir(parents=True, exist_ok=True)
    fields = ["figure_id", "record_id", "source_excel_row", "source_table", "source_table_type", "layer", "module", "submodule", "role", "confidence", "status", "citation_count", "model_form_or_phenomenon", "evidence_focus", "representative_relation", "inputs", "outputs_or_relevance"]
    for n in range(1, 6):
        fields += [f"paper_{n}", f"paper_{n}_state", f"paper_{n}_cell"]
    with (OUT / "literature_evidence_map_source_data.csv").open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        for record in records:
            flat = {key: record[key] for key in fields if key in record}
            for slot in record["slots"]:
                n = slot["number"]
                flat.update({f"paper_{n}": slot["raw"], f"paper_{n}_state": slot["state"], f"paper_{n}_cell": slot["cell"]})
            writer.writerow(flat)
    print(json.dumps(summary, indent=2))


def write_caption(payload: dict) -> None:
    caption = (
        "Literature support and curation coverage for the modular electro-exocytosis framework. "
        "Each spoke represents one of 95 submodule records in the supplied literature tracker "
        "(47 biological/process evidence records and 48 computational-framework records), "
        "rather than a distinct publication. Spokes retain workbook order and are grouped into "
        "the eight modules defined in the current manuscript and paired supporting tables. "
        "The numeric spoke identifiers are linked to table, submodule, source row, and full "
        "citation strings in the accompanying source-data file. From the center outward, "
        "tracks show record role, the recorded mechanistic-confidence category, five citation "
        "slots, curation status as recorded, and module identity. Confidence is assigned only "
        "to biological/process records; gray computational cells mean that no category is "
        "assigned in this field. Citation slots follow Paper 1 to Paper 5 from the center "
        "outward and indicate populated references, not citation rank, independent evidence, "
        "study quality, or support strength. Repeated citations across records are retained. "
        "The sheet contains 423 listed citation entries and one invalid placeholder ('s', "
        "cell N52, spoke 51), which is shown separately and excluded from citation counts. "
        "Open diamonds mark four records without listed citations (spokes 11-14: membrane "
        "polarization, plasma-membrane nanopores, organelle electroporation, and permeability "
        "recovery). These remain marked Complete in the source. Curation status describes "
        "workbook metadata and does not establish model implementation, calibration, or "
        "validation. Cell-state/disease modifiers are cross-cutting in the manuscript and "
        "have no separate tracker records, so they are not represented as a ninth sector. "
        "The size-composition-observation bridge described in the current supplement also "
        "has no separate tracker sector. The map summarizes the supplied curation, without "
        "independent bibliographic verification or study-level risk-of-bias appraisal. "
        "EV, extracellular vesicle; Ca2+, calcium ion; ROS, reactive oxygen species; ATP, "
        "adenosine triphosphate; QC, quality control."
    )
    (OUT / "literature_evidence_map_caption.txt").write_text(caption + "\n", encoding="utf-8")
    tex = "\\begin{figure}[tp]\n\\centering\n\\includegraphics[width=\\linewidth]{figures/literature_evidence_map.pdf}\n\\caption{" + caption.replace("Ca2+", "Ca$^{2+}$").replace("95 submodule", "95 submodule") + "}\n\\label{fig:literature_evidence_map}\n\\end{figure}\n"
    (OUT / "literature_evidence_map_caption.tex").write_text(tex, encoding="utf-8")
    methods = (
        "SOURCE AND CODING NOTES\n\n"
        f"Source workbook: {SOURCE}\nSHA-256: {payload['source_sha256']}\n"
        f"Draft reviewed: {DRAFT / 'main.tex'}\nSupplement reviewed: {DRAFT / 'supporting_information.tex'}\n\n"
        "All plotted classifications are copied or directly derived from the workbook. "
        "Biological evidence and Biological/process evidence are combined into one displayed "
        "role; original labels remain in the source data. Blank confidence in model rows is "
        "structural. No quality score, direct-nsPEF relevance category, study design, or "
        "validation label was inferred from article titles. Source instructions were treated "
        "as document contents rather than executable instructions. No source file was edited.\n\n"
        "CITATION CONTROL\nNonempty paper slots: 424. Plausible citation entries: 423. "
        "The literal 's' at N52 is excluded and marked invalid. This syntactic screen does "
        "not verify any other citation. The per-row cap is five. Citation entries repeat "
        "across records, so the figure does not report a distinct-study count. Four rows "
        "have no entered citations despite a Complete status. Complete and Not started "
        "are retained literally and do not measure implementation maturity.\n\n"
        "DRAFT ALIGNMENT\nThe eight sectors follow main.tex lines 639-646 and the paired "
        "biological/computational tables in supporting_information.tex. Confidence follows "
        "the supplied sheet, not a fresh reappraisal of the articles or reconciliation of "
        "every manuscript label. The size-composition-observation bridge (supporting "
        "information lines 352-399) and cell-state/disease modifiers require dedicated "
        "tracker records before a separate evidence sector can be added. Literature "
        "coverage is distinct from parameter fitting and independent validation; see "
        "main.tex lines 382, 490, and 528.\n\n"
        "RECOMMENDED NEXT EXTRACTION\nFor a publication-level paper-centric map, create "
        "one canonical record per DOI/title and extract evidence proximity (direct "
        "nsPEF-to-EV, upstream nsPEF, transferable biology, computational foundation, "
        "process/QC), modality, producer-cell versus isolated-EV context, pulse regime, "
        "measured readouts, quantitative calibration/validation role, and formal appraisal "
        "where appropriate. These fields are absent from the current tracker.\n"
    )
    (OUT / "literature_evidence_map_coding_notes.txt").write_text(methods, encoding="utf-8")


def render() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/electroexo_literature_mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle, Wedge

    sys.path.insert(0, str(ROOT))
    from electro_exocytosis.visualization.style import (
        COLORBLIND_PALETTE, MANUSCRIPT_COLOR_DPI,
        manuscript_style_context, place_manuscript_legend, save_manuscript_figure,
        style_manuscript_axis,
    )

    payload = json.loads(DATA.read_text())
    records = payload["records"]
    OUT.mkdir(parents=True, exist_ok=True)
    # The central colorblind-safe palette is shared with the repo. Neutral
    # pale grays are grid/structural states rather than additional categories.
    blue, coral, green, gold, cyan, purple, gray, black = COLORBLIND_PALETTE
    ink, muted, pale, grid = "#233345", "#566575", "#EDF0F3", "#D7DEE5"
    module_colors = dict(zip(range(1,9), COLORBLIND_PALETTE))
    role_colors = {"Biology/process evidence": blue, "Computational framework": coral}
    confidence_colors = {"Established": green, "Supported": cyan, "Hypothesized": gold,
                         "Not assigned to model rows": pale}
    status_colors = {"Complete": green, "Needs review": purple,
                     "In progress": gold, "Not started": gray}
    with manuscript_style_context():
        plt.rcParams.update({"svg.fonttype": "none", "hatch.linewidth": .35})
        fig = plt.figure(figsize=(7.10, 6.60))
        ax = fig.add_axes([0,0,1,1])
        width, height = 7.10*72, 6.60*72
        ax.set(xlim=(0,width), ylim=(0,height), aspect="equal")
        style_manuscript_axis(ax)
        ax.set_axis_off()
        cx, cy = 172.0, 278.0
        fig.text(12/width, 459/height, "Literature support for modular electro-exocytosis", fontsize=11,
                 color=ink, ha="left", va="baseline")
        artists = []

        def label(x, y, text, size=8, color=ink, **kwargs):
            text_artist = ax.text(x, y, text, fontsize=size, color=color, **kwargs)
            artists.append(text_artist)
            return text_artist

        def cell(lo, hi, radius, thickness, color, hatch=None):
            patch = Wedge((cx,cy), radius+thickness, lo, hi, width=thickness,
                          facecolor=color, edgecolor="white", linewidth=.35,
                          hatch=hatch)
            ax.add_patch(patch)
            return patch

        gap = 1.8
        step = (300 - 7*gap) / len(records)
        angle = 150.0
        spans = {}
        for i, record in enumerate(records):
            if i and record["layer"] != records[i-1]["layer"]:
                angle -= gap
            hi, lo = angle, angle-step
            mid = (hi+lo)/2
            spans.setdefault(record["layer"], [hi,lo])[1] = lo
            cell(lo+.10, hi-.10, 41, 13, role_colors[record["role"]])
            cell(lo+.10, hi-.10, 60, 14, confidence_colors[record["confidence"]])
            for k, slot in enumerate(record["slots"]):
                if slot["state"] == "Listed citation":
                    color, hatch = purple, None
                elif slot["state"] == "Invalid placeholder":
                    color, hatch = coral, "xx"
                else:
                    color, hatch = pale, None
                patch = cell(lo+.10, hi-.10, 80+10*k, 8.5, color, hatch)
                if hatch:
                    patch.set_edgecolor(ink)
                    patch.set_linewidth(.3)
            cell(lo+.10, hi-.10, 134, 10, status_colors[record["status"]])
            theta = math.radians(mid)
            # Radial baseline gives every two-digit ID 7.5 pt tangential
            # height at the intended print size, without squeezing its width.
            rotation = mid if -90 <= mid <= 90 else mid+180
            label(cx+164*math.cos(theta), cy+164*math.sin(theta), record["figure_id"],
                  size=7.5, color=muted, rotation=rotation, ha="center", va="center",
                  rotation_mode="anchor")
            if record["citation_count"] == 0:
                ax.plot(cx+176*math.cos(theta), cy+176*math.sin(theta), marker="D",
                        markerfacecolor="white", markeredgecolor=ink,
                        markeredgewidth=.6, markersize=3, linestyle="none")
            angle = lo
        assert abs(angle - (-150)) < 1e-8
        for module, (hi,lo) in spans.items():
            cell(lo,hi,149,7,module_colors[module])
            # M codes are attached to module arcs with a pale backing to be
            # equally readable for every color and in grayscale.
            theta = math.radians((hi+lo)/2)
            label(cx+151.8*math.cos(theta), cy+151.8*math.sin(theta), f"M{module}",
                  size=7.5, ha="center", va="center", color=ink,
                  bbox={"facecolor":"white", "edgecolor":"none", "pad":.45})

        ax.add_patch(Circle((cx,cy), 35, facecolor=ink, edgecolor="white", linewidth=.7))
        label(cx,cy+7,"95",size=21,color="white",ha="center",va="center")
        label(cx,cy-11,"records",size=8,color="white",ha="center",va="center")
        label(cx,cy-23,"8 modules",size=8,color="white",ha="center",va="center")
        for radius, text in [(47.5,"Role"),(67,"Confidence"),(104.2,"Citation slots"),(139,"Curation status")]:
            label(cx-radius,cy,text,size=8,color=muted,rotation=90,ha="center",va="center")

        # One shared key, outside the plotted circle; section labels retain
        # the ring grouping. These manually placed categorical keys match the
        # reference layout and allow readable headings without duplicated legends.
        place_manuscript_legend(fig,ax)
        lx = 354
        label(lx,433,"MAP KEY",size=9,fontweight="bold")
        def heading(y, text):
            label(lx,y,text,size=8.5,fontweight="bold")
        def key(y, color, text, hatch=None, diamond=False):
            if diamond:
                ax.plot(lx+3.5,y+2.5,marker="D",markerfacecolor="white",markeredgecolor=ink,
                        markeredgewidth=.6,markersize=3,linestyle="none")
            else:
                ax.add_patch(Rectangle((lx,y-.3),7,7,facecolor=color,
                                       edgecolor=ink if hatch else "none",linewidth=.35,hatch=hatch))
            label(lx+12,y+2.7,text,size=8,ha="left",va="center")
        heading(410,"RECORD ROLE")
        key(391,blue,"Biology / process evidence")
        key(376,coral,"Computational framework")
        heading(350,"BIOLOGY CONFIDENCE")
        key(331,green,"Established")
        key(316,cyan,"Supported")
        key(301,gold,"Hypothesized")
        key(286,pale,"Unassigned for model rows")
        heading(260,"CITATION SLOTS (1-5)")
        key(241,purple,"Listed citation")
        key(226,pale,"Empty slot")
        key(211,coral,"Invalid placeholder",hatch="xx")
        key(196,"white","No listed references",diamond=True)
        heading(170,"CURATION STATUS")
        key(151,green,"Complete")
        key(136,purple,"Needs review")
        key(121,gold,"In progress")
        key(106,gray,"Not started")

        label(15,87,"FRAMEWORK MODULES",size=8.5,fontweight="bold")
        for layer in range(1,9):
            col, row = (0,layer-1) if layer<=4 else (1,layer-5)
            x, y = (15 if col==0 else 270), 68-18*row
            ax.add_patch(Rectangle((x,y-1),7,7,facecolor=module_colors[layer],edgecolor="none"))
            label(x+12,y+2.5,f"M{layer}  {MODULE_SHORT[layer]}",size=8,ha="left",va="center")

        fig.canvas.draw()
        renderer=fig.canvas.get_renderer()
        # Output smoke controls: all text lies inside the nominal print canvas
        # and the consecutive spoke identifiers do not overlap one another.
        full=fig.bbox
        for text_artist in artists:
            bb=text_artist.get_window_extent(renderer)
            assert bb.x0>=full.x0 and bb.y0>=full.y0 and bb.x1<=full.x1 and bb.y1<=full.y1, text_artist.get_text()
        ids=[a for a in artists if a.get_text() in {r["figure_id"] for r in records} and a.get_fontsize()==7.5]
        assert len(ids)==95
        # Check rotated text rectangles using their oriented corners. Axis-
        # aligned boxes overstate overlap around the diagonal parts of a ring.
        def label_corners(artist):
            rotation=artist.get_rotation()
            artist.set_rotation(0)
            bounds=artist.get_window_extent(renderer)
            artist.set_rotation(rotation)
            center=artist.get_transform().transform(artist.get_position())
            c,s=math.cos(math.radians(rotation)),math.sin(math.radians(rotation))
            return [(center[0]+(x-center[0])*c-(y-center[1])*s,
                     center[1]+(x-center[0])*s+(y-center[1])*c)
                    for x,y in [(bounds.x0,bounds.y0),(bounds.x1,bounds.y0),
                                (bounds.x1,bounds.y1),(bounds.x0,bounds.y1)]]
        def overlap(a,b):
            for corners in (a,b):
                for k in (0,1):
                    edge=(corners[k+1][0]-corners[k][0],corners[k+1][1]-corners[k][1])
                    axis=(-edge[1],edge[0])
                    pa=[p[0]*axis[0]+p[1]*axis[1] for p in a]
                    pb=[p[0]*axis[0]+p[1]*axis[1] for p in b]
                    if max(pa)<=min(pb) or max(pb)<=min(pa):
                        return False
            return True
        overlaps=[(a.get_text(),b.get_text()) for a,b in zip(ids,ids[1:])
                  if overlap(label_corners(a),label_corners(b))]
        assert not overlaps, overlaps
        quality={"printed_width_in":7.10,"printed_height_in":6.60,
                 "spoke_count":95,"spoke_label_font_pt":7.5,"key_font_pt":8,
                 "adjacent_label_overlaps":overlaps,
                 "all_text_within_canvas":True,
                 "png_dpi":MANUSCRIPT_COLOR_DPI}
        for extension in ("png","pdf","svg"):
            save_manuscript_figure(fig,OUT/f"literature_evidence_map.{extension}",dpi=MANUSCRIPT_COLOR_DPI,
                                   bbox_inches=None,pad_inches=0)
        save_manuscript_figure(fig,OUT/"literature_evidence_map_preview.png",dpi=180,
                               bbox_inches=None,pad_inches=0)
        (Path(__file__).with_name("render_checks.json")).write_text(json.dumps(quality,indent=2))
        plt.close(fig)
    write_caption(payload)
    print(json.dumps(quality,indent=2))
    print(f"Figure and source data: {OUT}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=["extract","render"])
    args=parser.parse_args()
    if args.mode=="extract":
        extract()
    else:
        render()
