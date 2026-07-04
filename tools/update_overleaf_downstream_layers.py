#!/usr/bin/env python3
"""Patch the Overleaf manuscript with downstream Layer 6-8 implementation results."""

from __future__ import annotations

import shutil
from pathlib import Path


OVERLEAF = Path("/Users/aghorban/code/electro-exocytosis")
FIG_SRC = Path("/Users/aghorban/code/electroexo/results/downstream_engineering_comparison/downstream_engineering_panel.png")
FIG_DST = OVERLEAF / "figs" / "downstream_engineering_panel.png"


BIB_ENTRIES = r"""
@article{huotari2011endosome,
  author = {Huotari, Jatta and Helenius, Ari},
  title = {Endosome maturation},
  journal = {The EMBO Journal},
  year = {2011},
  volume = {30},
  number = {17},
  pages = {3481--3500},
  doi = {10.1038/emboj.2011.286}
}

@article{jeppesen2019reassessment,
  author = {Jeppesen, Dennis K. and Fenix, A. M. and Franklin, J. L. and Higginbotham, J. N. and Zhang, Q. and Zimmerman, L. J. and others},
  title = {Reassessment of exosome composition},
  journal = {Cell},
  year = {2019},
  volume = {177},
  number = {2},
  pages = {428--445.e18},
  doi = {10.1016/j.cell.2019.02.029}
}

@article{kugeratski2021quantitative,
  author = {Kugeratski, Fernanda G. and Hodge, Katie and Lilla, Stefania and McAndrews, Kathleen M. and Zhou, Xiao and Hwang, Robert F. and others},
  title = {Quantitative proteomics identifies the core proteome of exosomes with syntenin-1 as the highest abundant protein and a putative universal biomarker},
  journal = {Nature Cell Biology},
  year = {2021},
  volume = {23},
  number = {6},
  pages = {631--641},
  doi = {10.1038/s41556-021-00693-y}
}

@article{shurtleff2016ybox,
  author = {Shurtleff, Matthew J. and Temoche-Diaz, Margaret M. and Karfilis, Katherine V. and Ri, Shun and Schekman, Randy},
  title = {Y-box protein 1 is required to sort microRNAs into exosomes in cells and in a cell-free reaction},
  journal = {eLife},
  year = {2016},
  volume = {5},
  pages = {e19276},
  doi = {10.7554/eLife.19276}
}

@article{liu2021selective,
  author = {Liu, Xiao-Ming and Ma, Lu and Schekman, Randy},
  title = {Selective sorting of microRNAs into exosomes by phase-separated YBX1 condensates},
  journal = {eLife},
  year = {2021},
  volume = {10},
  pages = {e71982},
  doi = {10.7554/eLife.71982}
}

@article{santangelo2016syncrip,
  author = {Santangelo, Laura and Giurato, Giovanna and Cicchini, Claudia and Montaldo, Claudio and Mancone, Carmela and Tarallo, Roberta and others},
  title = {The RNA-binding protein SYNCRIP is a component of the hepatocyte exosomal machinery controlling microRNA sorting},
  journal = {Cell Reports},
  year = {2016},
  volume = {17},
  number = {3},
  pages = {799--808},
  doi = {10.1016/j.celrep.2016.09.031}
}

@article{skotland2023lipids,
  author = {Skotland, Tore and Llorente, Alicia and Sandvig, Kirsten},
  title = {Lipids in extracellular vesicles: what can be learned about membrane structure and function?},
  journal = {Cold Spring Harbor Perspectives in Biology},
  year = {2023},
  volume = {15},
  number = {8},
  pages = {a041415},
  doi = {10.1101/cshperspect.a041415}
}

@article{fuhrmann2015active,
  author = {Fuhrmann, Gregor and Serio, Andrea and Mazo, Maria and Nair, Rohit and Stevens, Molly M.},
  title = {Active loading into extracellular vesicles significantly improves the cellular uptake and photodynamic effect of porphyrins},
  journal = {Journal of Controlled Release},
  year = {2015},
  volume = {205},
  pages = {35--44},
  doi = {10.1016/j.jconrel.2014.11.029}
}

@article{didiot2018optimized,
  author = {Didiot, Marie-Cecile and Hall, Linda M. and Coles, Anne H. and Haraszti, Reka A. and Godinho, Bruno M. D. C. and Chase, Kris and others},
  title = {Optimized cholesterol-siRNA chemistry improves productive loading onto extracellular vesicles},
  journal = {Molecular Therapy},
  year = {2018},
  volume = {26},
  number = {6},
  pages = {1500--1509},
  doi = {10.1016/j.ymthe.2018.05.024}
}

@article{albeck2008modeling,
  author = {Albeck, John G. and Burke, John M. and Spencer, Sabrina L. and Lauffenburger, Douglas A. and Sorger, Peter K.},
  title = {Modeling a snap-action, variable-delay switch controlling extrinsic cell death},
  journal = {PLoS Biology},
  year = {2008},
  volume = {6},
  number = {12},
  pages = {e299},
  doi = {10.1371/journal.pbio.0060299}
}

@article{bentele2004mathematical,
  author = {Bentele, Marc and Lavrik, Inna and Ulrich, Marcel and Stosser, Stefanie and Heermann, Dieter W. and Kalthoff, Holger and others},
  title = {Mathematical modeling reveals threshold mechanism in CD95-induced apoptosis},
  journal = {The Journal of Cell Biology},
  year = {2004},
  volume = {166},
  number = {6},
  pages = {839--851},
  doi = {10.1083/jcb.200404158}
}

@article{cappe2023systematic,
  author = {Cappe, Benjamin and Vadi, Maija and Sack, Ethan and Wacheul, Laurent and Verstraeten, Bert and Dufour, Simon and others},
  title = {Systematic compositional analysis of exosomal extracellular vesicles produced by cells undergoing apoptosis, necroptosis and ferroptosis},
  journal = {Journal of Extracellular Vesicles},
  year = {2023},
  volume = {12},
  number = {10},
  pages = {e12365},
  doi = {10.1002/jev2.12365}
}

@article{webber2013pure,
  author = {Webber, James and Clayton, Aled},
  title = {How pure are your vesicles?},
  journal = {Journal of Extracellular Vesicles},
  year = {2013},
  volume = {2},
  number = {1},
  pages = {19861},
  doi = {10.3402/jev.v2i0.19861}
}

@article{brennan2020comparison,
  author = {Brennan, K. and Martin, K. and FitzGerald, S. P. and O'Sullivan, J. and Wu, Y. and Blanco, A. and others},
  title = {A comparison of methods for the isolation and separation of extracellular vesicles from protein and lipid particles in human serum},
  journal = {Scientific Reports},
  year = {2020},
  volume = {10},
  pages = {1039},
  doi = {10.1038/s41598-020-57497-7}
}

@article{paganini2019scalable,
  author = {Paganini, Cristina and Capasso Palmiero, Ugo and Pocsfalvi, Gabriella and Touzet, Nicolas and Bongiovanni, Antonella and Arosio, Paolo},
  title = {Scalable production and isolation of extracellular vesicles: available sources and lessons from current industrial bioprocesses},
  journal = {Biotechnology Journal},
  year = {2019},
  volume = {14},
  number = {10},
  pages = {1800528},
  doi = {10.1002/biot.201800528}
}

@article{watson2018scalable,
  author = {Watson, Douglas C. and Yung, Benjamin C. and Bergamaschi, Cristina and Chowdhury, Biswajit and Bear, Joseph and Stellas, Dimitrios and others},
  title = {Scalable, cGMP-compatible purification of extracellular vesicles carrying bioactive human heterodimeric IL-15/lactadherin complexes},
  journal = {Journal of Extracellular Vesicles},
  year = {2018},
  volume = {7},
  number = {1},
  pages = {1442088},
  doi = {10.1080/20013078.2018.1442088}
}

@article{mendt2018generation,
  author = {Mendt, Mayela and Kamerkar, Shilpa and Sugimoto, Hiroshi and McAndrews, Kathleen M. and Wu, Chien-Chia and Gagea, Mihai and others},
  title = {Generation and testing of clinical-grade exosomes for pancreatic cancer},
  journal = {JCI Insight},
  year = {2018},
  volume = {3},
  number = {8},
  pages = {e99263},
  doi = {10.1172/jci.insight.99263}
}

@article{thakur2024global,
  author = {Thakur, Amit and Rai, Devendra},
  title = {Global requirements for manufacturing and validation of clinical grade extracellular vesicles},
  journal = {The Journal of Liquid Biopsy},
  year = {2024},
  volume = {6},
  pages = {100278},
  doi = {10.1016/j.jlb.2024.100278}
}
"""


PRELIM_INSERT_AFTER = """The integrated storyboard provides a compact narrative view of the framework. In the mild reversible window, the model predicts threshold-level charging, modest temperature rise, limited Ca$^{2+}$ amplification, and predominantly repair-compatible remodeling with low apoptotic output. In the productive secretory window, larger Ca$^{2+}$ loading and stronger repair support coexist with a more favorable secretory routing bias, increasing small-EV output without a catastrophic loss of viability. In the injury-dominant window, the same layered structure yields a different interpretation: stronger charging and thermal load amplify Ca$^{2+}$ and ROS, repair remains incomplete, and apoptotic vesiculation becomes a dominant endpoint. This row-wise transition is encouraging because it shows that the current framework already produces coherent scenario-level behavior while remaining modular enough for later biological calibration.
"""

PRELIM_INSERT = r"""
\begin{figure*}[t]
\centering
\includegraphics[width=0.98\textwidth]{figs/downstream_engineering_panel.png}
\caption{\textbf{Preliminary Layer 6--8 downstream engineering behavior.} The downstream comparison extends the preliminary examples from vesicle release into product interpretation. A productive cargo window increases protein, RNA, lipid, and antigen cargo scores through secretory and ESCRT-linked terms; a direct-loading mode adds a permeability-driven cargo term while retaining acceptable viability; and an injury-dominant condition may still produce particles but fails the quality gate because apoptosis, necrosis, and contaminant fractions reduce bona fide EV purity and process objective.}
\abbrevnote{EV, extracellular vesicle; RNA, ribonucleic acid; ESCRT, endosomal sorting complex required for transport.}
\label{fig:prelim_downstream_engineering}
\end{figure*}

The Layer 6--8 comparison illustrates why the implementation now carries cargo, quality, and manufacturing states beyond total particle release. Productive and direct-loading scenarios can yield similar process objectives through different mechanisms--endogenous cargo sorting versus permeability-assisted loading--whereas the injury-dominant scenario shows that a high particle output is not automatically a useful EV product. This behavior follows the report's computational tables: subtype-weighted cargo and stress sorting feed potency, damage and contaminant mixtures determine the quality gate, and recovery, purity, batch consistency, and potency jointly determine the manufacturing objective.
"""

LAYER6_AFTER = r"""\abbrevnote{EV, extracellular vesicle; ESCRT, endosomal sorting complex required for transport; ALIX, ALG-2-interacting protein X; sEV, small extracellular vesicle; m/lEV, medium/large extracellular vesicle; AB, apoptotic body.}
\end{table*}
"""

LAYER6_TEXT = r"""

The implemented Layer 6 module now follows the table structure explicitly. It computes a subtype-weighted cargo state from small-EV, medium/large-EV, and apoptotic-body outputs; separates protein, RNA, lipid, antigen, and direct-loading terms; and maps the resulting cargo vector to a saturating potency score. This keeps the implementation aligned with evidence that EV subtypes carry distinct protein compositions \cite{jeppesen2019reassessment,kugeratski2021quantitative}, that RNA sorting can be mediated by RNA-binding proteins such as YBX1 and SYNCRIP \cite{shurtleff2016ybox,liu2021selective,santangelo2016syncrip}, that lipid composition and ceramide biology shape EV state \cite{skotland2023lipids}, and that direct EV loading can be represented separately from endogenous cargo sorting \cite{fuhrmann2015active,didiot2018optimized}.
"""

LAYER7_AFTER = r"""\abbrevnote{EV, extracellular vesicle; AB, apoptotic body; ATP, adenosine triphosphate; ROS, reactive oxygen species; QC, quality control.}
\end{table*}
"""

LAYER7_TEXT = r"""

The implemented Layer 7 module expands the previous binary quality proxy into a mixture-aware gate. Damage is mapped to stressed-viable, apoptotic, and necrotic fractions using threshold-like Hill terms inspired by systems models of cell-death switching \cite{albeck2008modeling,bentele2004mathematical}. The measured particle pool is then decomposed into bona fide EVs, apoptotic bodies, debris, and aggregate-like contaminants, which makes the purity score sensitive to both injury state and non-EV co-isolates. This structure matches EV-composition and isolation cautions from apoptotic/necroptotic EV studies and MISEV-style purity guidance \cite{cappe2023systematic,webber2013pure,brennan2020comparison,welsh2024minimal}.
"""

LAYER8_AFTER = r"""\abbrevnote{EV, extracellular vesicle; nsPEF, nanosecond pulsed electric field; NTA, nanoparticle tracking analysis; ELISA, enzyme-linked immunosorbent assay; QC, quality control; TEM, transmission electron microscopy.}
\end{table*}
"""

LAYER8_TEXT = r"""

The implemented Layer 8 module converts biological output into process-facing quantities: cell-normalized yield, harvest-rate-equivalent accumulation, isolation recovery, process purity, batch-adjusted yield, and a constrained optimization objective. This keeps the code aligned with EV manufacturing literature emphasizing scalable source material, recovery and purity tradeoffs, clinical-grade purification, and validation requirements \cite{paganini2019scalable,watson2018scalable,mendt2018generation,thakur2024global,welsh2024minimal}. The objective is intentionally transparent rather than prescriptive: potency, yield, purity, and viability weights can be changed for different process-development priorities.
"""


def append_missing_bib_entries(bib_path: Path) -> None:
    text = bib_path.read_text(encoding="utf-8")
    additions = []
    for block in BIB_ENTRIES.strip().split("\n\n@"):
        entry = block if block.startswith("@") else "@" + block
        key = entry.split("{", 1)[1].split(",", 1)[0]
        if f"{{{key}," not in text:
            additions.append(entry)
    if additions:
        bib_path.write_text(text.rstrip() + "\n\n" + "\n\n".join(additions) + "\n", encoding="utf-8")


def insert_once(text: str, anchor: str, insertion: str, marker: str) -> str:
    if marker in text:
        return text
    if anchor not in text:
        raise RuntimeError(f"anchor not found for marker {marker}")
    return text.replace(anchor, anchor + insertion, 1)


def patch_main(main_path: Path) -> None:
    text = main_path.read_text(encoding="utf-8")
    text = insert_once(text, PRELIM_INSERT_AFTER, PRELIM_INSERT, "fig:prelim_downstream_engineering")
    text = insert_once(text, LAYER6_AFTER, LAYER6_TEXT, "The implemented Layer 6 module now follows")
    text = insert_once(text, LAYER7_AFTER, LAYER7_TEXT, "The implemented Layer 7 module expands")
    text = insert_once(text, LAYER8_AFTER, LAYER8_TEXT, "The implemented Layer 8 module converts")
    main_path.write_text(text, encoding="utf-8")


def main() -> int:
    FIG_DST.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(FIG_SRC, FIG_DST)
    append_missing_bib_entries(OVERLEAF / "source.bib")
    patch_main(OVERLEAF / "main.tex")
    print(f"copied {FIG_DST}")
    print("patched main.tex and source.bib")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
