#!/usr/bin/env python3
"""Add row-level citations to Overleaf appendix tables from the evidence tracker."""

from __future__ import annotations

from pathlib import Path


OVERLEAF = Path("/Users/aghorban/code/electro-exocytosis")
MAIN = OVERLEAF / "main.tex"
BIB = OVERLEAF / "source.bib"


BIB_ENTRIES = r"""
@article{singh2025electroporation,
  author = {Singh, M. and Mazaheri-Tehrani, G. and Martin-Fabiani, I. and Davies, O. G.},
  title = {Electroporation induced changes in extracellular vesicle profile},
  journal = {Drug Delivery},
  year = {2025},
  volume = {32},
  number = {1},
  pages = {2562224},
  doi = {10.1080/10717544.2025.2562224}
}

@article{hood2014maximizing,
  author = {Hood, Joshua L. and Scott, Mark J. and Wickline, Samuel A.},
  title = {Maximizing exosome colloidal stability following electroporation},
  journal = {Analytical Biochemistry},
  year = {2014},
  volume = {448},
  pages = {41--49},
  doi = {10.1016/j.ab.2013.12.001}
}

@article{kooijmans2013electroporation,
  author = {Kooijmans, Sander A. A. and Stremersch, Sam and Braeckmans, Kevin and de Smedt, Stefaan C. and Hendrix, An and Wood, Matthew J. A. and Schiffelers, Raymond M. and Raemdonck, Koen and Vader, Pieter},
  title = {Electroporation-induced siRNA precipitation obscures the efficiency of siRNA loading into extracellular vesicles},
  journal = {Journal of Controlled Release},
  year = {2013},
  volume = {172},
  number = {1},
  pages = {229--238},
  doi = {10.1016/j.jconrel.2013.08.014}
}

@article{szlasa2024mage,
  author = {Szlasa, Wojciech and Sauer, Natalia and Baczynska, Dominika and Zietek, Michal and Haczkiewicz-Lesniak, Katarzyna and Karpinski, Piotr and Fleszar, Marcin and Fortuna, Filip and Kulus, Marta J. and Piotrowska, Aleksandra and Kmiecik, Aleksandra and Baranska, Agata and Michel, Olivier and Novickij, Vitalij and Tarek, Mounir and Kasperkiewicz, Paulina and Dziegiel, Piotr and Podhorska-Okolow, Marzena and Saczko, Jolanta and Kulbacka, Julita},
  title = {Pulsed electric field induces exocytosis and overexpression of MAGE antigens in melanoma},
  journal = {Scientific Reports},
  year = {2024},
  volume = {14},
  pages = {12546},
  doi = {10.1038/s41598-024-63181-x}
}

@article{kumar2022tumour,
  author = {Kumar, Pankaj and Boyne, Connor and Brown, Samantha and Qureshi, Arif and Thorpe, Paul and Synowsky, Sander A. and others},
  title = {Tumour-associated antigenic peptides are present in the HLA class I ligandome of cancer cell line derived extracellular vesicles},
  journal = {Immunology},
  year = {2022},
  volume = {166},
  number = {2},
  pages = {249--264},
  doi = {10.1111/imm.13471}
}

@article{hagey2023cellular,
  author = {Hagey, Daniel W. and Ojansivu, Miina and Bostancioglu, Burcu R. and Saher, Osama and Bost, Jean-Philippe and Gustafsson, Maria O. and Gramignoli, Roberto and Svahn, Helene Andersson and Gupta, Divya and Stevens, Molly M. and Gorgens, Andre and El Andaloussi, Samir},
  title = {The cellular response to extracellular vesicles is dependent on their cell source and dose},
  journal = {Science Advances},
  year = {2023},
  volume = {9},
  number = {35},
  pages = {eadh1168},
  doi = {10.1126/sciadv.adh1168}
}

@article{tkach2017qualitative,
  author = {Tkach, Meriem and Kowal, Joanna and Zucchetti, Anna E. and Enserink, Lara and Jouve, Mickael and Lankar, Daniel and Saitakis, Michael and Martin-Jaular, Lorena and Thery, Clotilde},
  title = {Qualitative differences in T-cell activation by dendritic cell-derived extracellular vesicle subtypes},
  journal = {The EMBO Journal},
  year = {2017},
  volume = {36},
  number = {20},
  pages = {3012--3028},
  doi = {10.15252/embj.201696003}
}
"""


REPLACEMENTS = {
    # A2
    r"Pulse input & External forcing function & $E(t)=E_0u(t;t_p,N,f)$ & amplitude, width, number, repetition & $E(t)$ \\": r"Pulse input & External forcing function \cite{beebe2003nanosecond,vernier2003calcium,pakhomov2007membrane,semenov2013primary} & $E(t)=E_0u(t;t_p,N,f)$ & amplitude, width, number, repetition & $E(t)$ \\",
    r"Dose metric & Selectable absorbed-energy model & $u=\sigma E_{\mathrm{peak}}^2t_pN\chi_w\langle g_E^2\rangle$ & pulse train, conductivity, waveform, geometry & energy density, dose index \\": r"Dose metric & Selectable absorbed-energy model \cite{orlacchio2023effects,gudvangen2022electroporation} & $u=\sigma E_{\mathrm{peak}}^2t_pN\chi_w\langle g_E^2\rangle$ & pulse train, conductivity, waveform, geometry & energy density, dose index \\",
    r"Field distribution & Laplace/quasi-static field model & $\nabla\cdot(\sigma\nabla\phi)=0$, $E=-\nabla\phi$ & chamber geometry, conductivity & local $E(\mathbf{x},t)$ \\": r"Field distribution & Laplace/quasi-static field model \cite{joshi2000electroporation,kotnik2019membrane,baker2024numerical,leveque2017measurement} & $\nabla\cdot(\sigma\nabla\phi)=0$, $E=-\nabla\phi$ & chamber geometry, conductivity & local $E(\mathbf{x},t)$ \\",
    r"Thermal constraint & Adiabatic or lumped heat-loss model & $\Delta T_{\mathrm{end}}=F_{\mathrm{ret}}\Delta T_{\mathrm{ad}}$ & field, medium, cooling, pulse train & $T(t)$, $\Delta T_{\mathrm{end}}$ \\": r"Thermal constraint & Adiabatic or lumped heat-loss model \cite{orlacchio2023effects,gudvangen2022electroporation,kanduser2008temperature,song2011synergistic,song2017thermal} & $\Delta T_{\mathrm{end}}=F_{\mathrm{ret}}\Delta T_{\mathrm{ad}}$ & field, medium, cooling, pulse train & $T(t)$, $\Delta T_{\mathrm{end}}$ \\",
    r"Flow-through exposure & Residence-time model & $D_E(\mathbf{x})=\int_{\tau} E(\mathbf{x}(t))^2dt$ & flow rate, chamber field & exposure distribution \\": r"Flow-through exposure & Residence-time model \cite{liu2026microfluidic,garcia2016microfluidic,bonakdar2017microfluidic,li2015electroporation} & $D_E(\mathbf{x})=\int_{\tau} E(\mathbf{x}(t))^2dt$ & flow rate, chamber field & exposure distribution \\",
    # A3/A4
    r"Isolated EV permeabilization & Direct nsPEF effect on EV membranes & Vesicle integrity, cargo loading, aggregation, marker retention & Defines second mode for post-isolation EV engineering & Hypothesized; not implemented in the current cell-based layer \\": r"Isolated EV permeabilization & Direct nsPEF effect on EV membranes & Vesicle integrity, cargo loading, aggregation, marker retention & Defines second mode for post-isolation EV engineering & Hypothesized; not implemented in the current cell-based layer \cite{singh2025electroporation,hood2014maximizing,kooijmans2013electroporation,fuhrmann2015active,didiot2018optimized} \\",
    r"Membrane charging & Reduced Schwan/circuit model & $\Delta V_m=f_sE_{\mathrm{eff}}a(1-e^{-t_p/\tau_m})$ & $E(t)$, geometry, cell radius, membrane state & $\Delta V_m$ \\": r"Membrane charging & Reduced Schwan/circuit model \cite{kotnik2000analytical,joshi2000electroporation,kotnik2019membrane} & $\Delta V_m=f_sE_{\mathrm{eff}}a(1-e^{-t_p/\tau_m})$ & $E(t)$, geometry, cell radius, membrane state & $\Delta V_m$ \\",
    r"Organelle charging & Reduced multi-shell surrogate & $\Delta V_o=\alpha_o\Delta V_m$ & organelle factors, plasma-membrane voltage & $\Delta V_o$ \\": r"Organelle charging & Reduced multi-shell surrogate \cite{semenov2013primary,beebe2012transient,baker2024numerical} & $\Delta V_o=\alpha_o\Delta V_m$ & organelle factors, plasma-membrane voltage & $\Delta V_o$ \\",
    r"Isolated EV charging & Vesicle-scale circuit model & $C_{EV}\frac{d\Delta V_{EV}}{dt}+G_{EV}\Delta V_{EV}=f(E,r_{EV})$ & EV radius, medium conductivity & $P_{EV}(t)$, integrity state \\": r"Isolated EV charging & Vesicle-scale circuit model \cite{singh2025electroporation,hood2014maximizing,kooijmans2013electroporation} & $C_{EV}\frac{d\Delta V_{EV}}{dt}+G_{EV}\Delta V_{EV}=f(E,r_{EV})$ & EV radius, medium conductivity & $P_{EV}(t)$, integrity state \\",
    # A8
    r"Generic Ca$^{2+}$ sensor & Hill activation & $A_x=\frac{C_{sub}^n}{K_x^n+C_{sub}^n}$ & local Ca$^{2+}$, sensitivity & activated fraction \\": r"Generic Ca$^{2+}$ sensor & Hill activation \cite{naraghi1997linearized,shannon2004mathematical,scheffer2014mechanism} & $A_x=\frac{C_{sub}^n}{K_x^n+C_{sub}^n}$ & local Ca$^{2+}$, sensitivity & activated fraction \\",
    r"Calpain/cytoskeleton & Protease activation and substrate loss & $\frac{dF_{actin}}{dt}=-k_{calp}A_{calp}F_{actin}+k_{poly}$ & Ca$^{2+}$, calpain state & cortical remodeling \\": r"Calpain/cytoskeleton & Protease activation and substrate loss \cite{campbell2012structure,roberts2020ca,williams2025calpains} & $\frac{dF_{actin}}{dt}=-k_{calp}A_{calp}F_{actin}+k_{poly}$ & Ca$^{2+}$, calpain state & cortical remodeling \\",
    # A10
    r"MVB pool & Compartmental pool balance & $\frac{dM}{dt}=k_{mat}E_{endo}-k_{dock}M-k_{deg}M$ & endosome pool, maturation & MVB number \\": r"MVB pool & Compartmental pool balance \cite{huotari2011endosome,ott2025coordination,borchers2023regulatory,solinger2025escrting} & $\frac{dM}{dt}=k_{mat}E_{endo}-k_{dock}M-k_{deg}M$ & endosome pool, maturation & MVB number \\",
    r"ILV loading & Production-rate law & $\frac{dI}{dt}=k_{ILV}f_{ESCRT}M-k_{release}I$ & ESCRT/ceramide state, MVB pool & ILV load \\": r"ILV loading & Production-rate law \cite{takahashi2015hrs,wollert2010molecular,pfitzner2020escrt,larios2020alix,choezom2022nsmase2,crivelli2022ceramideTransfer} & $\frac{dI}{dt}=k_{ILV}f_{ESCRT}M-k_{release}I$ & ESCRT/ceramide state, MVB pool & ILV load \\",
    r"MVB docking & Trafficking/docking ODE & $\frac{dM_d}{dt}=k_{dock}f_{Rab}M-k_{fus}M_d$ & Rab activity, actin barrier & docked MVB pool \\": r"MVB docking & Trafficking/docking ODE \cite{hsu2010rab35,messenger2018munc13,izumi2021rab27,tang2026rab27,marie2023accessory} & $\frac{dM_d}{dt}=k_{dock}f_{Rab}M-k_{fus}M_d$ & Rab activity, actin barrier & docked MVB pool \\",
    r"Small-EV release & Ca$^{2+}$-sensitive fusion rate & $R_{sEV}=k_{fus}f_{Ca}(C_{sub})M_d$ & docked MVBs, Ca$^{2+}$, SNARE & small-EV release rate \\": r"Small-EV release & Ca$^{2+}$-sensitive fusion rate \cite{liu2023snare,matsui2023vamp5,woo2017munc13,courtney2023synaptotagmin7} & $R_{sEV}=k_{fus}f_{Ca}(C_{sub})M_d$ & docked MVBs, Ca$^{2+}$, SNARE & small-EV release rate \\",
    r"Budding-site formation & Budding pool ODE & $\frac{dB_{m/lEV}}{dt}=k_{bud}f_{PS}f_{tension}-k_{scission}B_{m/lEV}$ & PS, curvature, tension & budding precursor pool \\": r"Budding-site formation & Budding pool ODE \cite{pang2025arf,dai2019rhoa,menck2017nsmaseMicrovesicles} & $\frac{dB_{m/lEV}}{dt}=k_{bud}f_{PS}f_{tension}-k_{scission}B_{m/lEV}$ & PS, curvature, tension & budding precursor pool \\",
    r"Medium/large EV release & Scission-rate law & $R_{m/lEV}=k_{scission}f_{actomyosin}B_{m/lEV}$ & calpain, ROCK, actomyosin & medium/large EV rate \\": r"Medium/large EV release & Scission-rate law \cite{pang2025arf,dai2019rhoa,menck2017nsmaseMicrovesicles,williams2025calpains} & $R_{m/lEV}=k_{scission}f_{actomyosin}B_{m/lEV}$ & calpain, ROCK, actomyosin & medium/large EV rate \\",
    r"Apoptotic vesicle release & Injury-state weighted rate & $R_{AB}=k_{AB}F_{apop}$ & apoptotic fraction, blebbing rate & apoptotic-body release \\": r"Apoptotic vesicle release & Injury-state weighted rate \cite{teixeira2020rock1,ozdemir2021rock1,sakamaki2025mprip,cappe2023systematic} & $R_{AB}=k_{AB}F_{apop}$ & apoptotic fraction, blebbing rate & apoptotic-body release \\",
    # A11/A12
    r"Protein sorting & Ubiquitin/ESCRT/tetraspanin-dependent protein enrichment & ALIX, TSG101, CD63/CD81/CD9, antigen enrichment & Determines therapeutic or immunogenic protein cargo & Established \\": r"Protein sorting & Ubiquitin/ESCRT/tetraspanin-dependent protein enrichment & ALIX, TSG101, CD63/CD81/CD9, antigen enrichment & Determines therapeutic or immunogenic protein cargo & Established \cite{jeppesen2019reassessment,kugeratski2021quantitative,takahashi2015hrs,larios2020alix,rubinstein2025tetraspanins} \\",
    r"RNA sorting & RNA-binding protein-mediated miRNA/mRNA loading & hnRNPA2B1, SYNCRIP, YBX1, motif-based sorting & Enables nucleic-acid cargo modulation & Supported \\": r"RNA sorting & RNA-binding protein-mediated miRNA/mRNA loading & hnRNPA2B1, SYNCRIP, YBX1, motif-based sorting & Enables nucleic-acid cargo modulation & Supported \cite{shurtleff2016ybox,liu2021selective,santangelo2016syncrip} \\",
    r"Lipid sorting & Ceramide, cholesterol, sphingolipid, PE/PS enrichment & Lipidomics, nSMase2, tetraspanin domains & Controls EV stability, curvature, uptake, immunologic function & Established \\": r"Lipid sorting & Ceramide, cholesterol, sphingolipid, PE/PS enrichment & Lipidomics, nSMase2, tetraspanin domains & Controls EV stability, curvature, uptake, immunologic function & Established \cite{skotland2023lipids,choezom2022nsmase2,crivelli2022ceramideTransfer,rubinstein2025tetraspanins} \\",
    r"Stress-responsive cargo & ROS, hypoxia, inflammatory, apoptotic cargo changes & proteomics/lipidomics/RNA-seq after nsPEF & Distinguishes programmed cargo modulation from debris & Supported \\": r"Stress-responsive cargo & ROS, hypoxia, inflammatory, apoptotic cargo changes & proteomics/lipidomics/RNA-seq after nsPEF & Distinguishes programmed cargo modulation from debris & Supported \cite{pakhomova2012oxidative,cappe2023systematic,muratori2017activation} \\",
    r"Antigen enrichment & Tumor antigen expression and EV loading & MAGE and other antigen studies & Supports vaccine-oriented EV engineering & Supported \\": r"Antigen enrichment & Tumor antigen expression and EV loading & MAGE and other antigen studies & Supports vaccine-oriented EV engineering & Supported \cite{szlasa2024mage,kumar2022tumour} \\",
    r"Functional potency & Recipient-cell response & uptake, cargo delivery, immune activation, repair assays & Connects cargo state to biological activity & Hypothesized \\": r"Functional potency & Recipient-cell response & uptake, cargo delivery, immune activation, repair assays & Connects cargo state to biological activity & Hypothesized \cite{hagey2023cellular,tkach2017qualitative,fuhrmann2015active} \\",
    r"Subtype-weighted cargo & Weighted output model & $Cargo_j=w_{j,sEV}R_{sEV}+w_{j,m/lEV}R_{m/lEV}+w_{j,AB}R_{AB}$ & subtype rates, cargo weights & cargo abundance \\": r"Subtype-weighted cargo & Weighted output model \cite{jeppesen2019reassessment,kugeratski2021quantitative,cappe2023systematic} & $Cargo_j=w_{j,sEV}R_{sEV}+w_{j,m/lEV}R_{m/lEV}+w_{j,AB}R_{AB}$ & subtype rates, cargo weights & cargo abundance \\",
    r"Stress-modulated sorting & Rate law & $\frac{dCargo_j}{dt}=k_{sort,j}f_{stress}f_{Ca}f_{ESCRT}-k_{loss,j}Cargo_j$ & Ca$^{2+}$, ROS, ESCRT, stress & cargo state \\": r"Stress-modulated sorting & Rate law \cite{shurtleff2016ybox,liu2021selective,santangelo2016syncrip,pakhomova2012oxidative} & $\frac{dCargo_j}{dt}=k_{sort,j}f_{stress}f_{Ca}f_{ESCRT}-k_{loss,j}Cargo_j$ & Ca$^{2+}$, ROS, ESCRT, stress & cargo state \\",
    r"Antigen enrichment & Transcription/sorting coupling & $\frac{dA_g}{dt}=k_{expr}f_{stress}-k_{deg}A_g-k_{sort}A_g$ & stress, antigen expression & antigen cargo \\": r"Antigen enrichment & Transcription/sorting coupling \cite{szlasa2024mage,kumar2022tumour} & $\frac{dA_g}{dt}=k_{expr}f_{stress}-k_{deg}A_g-k_{sort}A_g$ & stress, antigen expression & antigen cargo \\",
    r"Lipid composition & Composition-state ODE & $\frac{dL_j}{dt}=k_{syn,j}f_{stress}-k_{sort,j}L_j$ & lipid metabolism, EV subtype & lipid cargo class \\": r"Lipid composition & Composition-state ODE \cite{skotland2023lipids,choezom2022nsmase2,crivelli2022ceramideTransfer} & $\frac{dL_j}{dt}=k_{syn,j}f_{stress}-k_{sort,j}L_j$ & lipid metabolism, EV subtype & lipid cargo class \\",
    r"Potency mapping & Empirical or mechanistic response model & $Potency=F(Cargo,Subtype,Dose)$ & cargo vector, recipient assay & functional potency \\": r"Potency mapping & Empirical or mechanistic response model \cite{hagey2023cellular,tkach2017qualitative,fuhrmann2015active} & $Potency=F(Cargo,Subtype,Dose)$ & cargo vector, recipient assay & functional potency \\",
    r"Direct EV loading & Permeability-diffusion model & $\frac{dC_{load}}{dt}=P_{EV}(t)A_{EV}(C_{ext}-C_{EV})-k_{leak}C_{EV}$ & EV permeability, cargo gradient & loaded cargo \\": r"Direct EV loading & Permeability-diffusion model \cite{kooijmans2013electroporation,hood2014maximizing,fuhrmann2015active,didiot2018optimized} & $\frac{dC_{load}}{dt}=P_{EV}(t)A_{EV}(C_{ext}-C_{EV})-k_{leak}C_{EV}$ & EV permeability, cargo gradient & loaded cargo \\",
    # A13/A14
    r"Viability preservation & Reversible response versus cell death & viability, membrane recovery, clonogenic/metabolic assays & Defines acceptable EV-production window & Established \\": r"Viability preservation & Reversible response versus cell death & viability, membrane recovery, clonogenic/metabolic assays & Defines acceptable EV-production window & Established \cite{pakhomova2013two,beebe2012transient,bhattacharya2022calcium} \\",
    r"Mitochondrial injury & Depolarization, mPTP, cytochrome c & $\Delta\Psi_m$, cytochrome c, ROS & Separates stress signaling from organellar contamination & Supported \\": r"Mitochondrial injury & Depolarization, mPTP, cytochrome c & $\Delta\Psi_m$, cytochrome c, ROS & Separates stress signaling from organellar contamination & Supported \cite{beebe2012transient,napotnik2012mitochondrial,pakhomova2012oxidative} \\",
    r"Apoptosis & Caspase activation and apoptotic blebbing & caspase, annexin V/PI, ROCK1 blebbing & Determines apoptotic-body contribution & Established \\": r"Apoptosis & Caspase activation and apoptotic blebbing & caspase, annexin V/PI, ROCK1 blebbing & Determines apoptotic-body contribution & Established \cite{pakhomova2013two,teixeira2020rock1,cappe2023systematic} \\",
    r"Necrosis/lysis & Irreversible membrane failure & LDH release, PI uptake, debris microscopy & Indicates unacceptable product contamination & Established \\": r"Necrosis/lysis & Irreversible membrane failure & LDH release, PI uptake, debris microscopy & Indicates unacceptable product contamination & Established \cite{pakhomova2013two,cappe2023systematic} \\",
    r"Non-EV contaminants & Protein aggregates, organelle fragments, apoptotic bodies & calnexin, cytochrome c, albumin/lipoprotein markers & Defines purity and interpretability of EV readouts & Established \\": r"Non-EV contaminants & Protein aggregates, organelle fragments, apoptotic bodies & calnexin, cytochrome c, albumin/lipoprotein markers & Defines purity and interpretability of EV readouts & Established \cite{webber2013pure,brennan2020comparison,cappe2023systematic,welsh2024minimal} \\",
    r"Quality-control marker panel & Positive and negative EV markers & CD63/CD81/CD9, TSG101/ALIX/syntenin, negative markers & Connects model to MISEV-style QC & Established \\": r"Quality-control marker panel & Positive and negative EV markers & CD63/CD81/CD9, TSG101/ALIX/syntenin, negative markers & Connects model to MISEV-style QC & Established \cite{welsh2024minimal,jeppesen2019reassessment,kugeratski2021quantitative} \\",
    r"Damage accumulation & Damage-state ODE & $\frac{dD}{dt}=k_{inj}G_m+k_{ROS}ROS-k_{repair}A_{repair}D$ & conductance, ROS, repair & damage burden \\": r"Damage accumulation & Damage-state ODE \cite{pakhomova2012oxidative,bhattacharya2022calcium,albeck2008modeling} & $\frac{dD}{dt}=k_{inj}G_m+k_{ROS}ROS-k_{repair}A_{repair}D$ & conductance, ROS, repair & damage burden \\",
    r"Fate transition & Hybrid state rule & $F_{apop}=H(D-D^*)$ or $\frac{dF_{apop}}{dt}=k_{apop}H(D-D^*)$ & damage, ATP, mitochondria & viable/stressed/apoptotic fractions \\": r"Fate transition & Hybrid state rule \cite{albeck2008modeling,bentele2004mathematical,pakhomova2013two} & $F_{apop}=H(D-D^*)$ or $\frac{dF_{apop}}{dt}=k_{apop}H(D-D^*)$ & damage, ATP, mitochondria & viable/stressed/apoptotic fractions \\",
    r"Necrotic transition & Threshold rule & $F_{nec}=H(D-D_{nec}^*)H(A_{ATP}^*-ATP)$ & damage, ATP, membrane failure & necrotic fraction \\": r"Necrotic transition & Threshold rule \cite{pakhomova2013two,cappe2023systematic} & $F_{nec}=H(D-D_{nec}^*)H(A_{ATP}^*-ATP)$ & damage, ATP, membrane failure & necrotic fraction \\",
    r"Measured particles & Mixture model & $N_{meas}=N_{EV}+N_{AB}+N_{debris}+N_{agg}$ & release rates, injury state & measured particle count \\": r"Measured particles & Mixture model \cite{cappe2023systematic,webber2013pure,brennan2020comparison} & $N_{meas}=N_{EV}+N_{AB}+N_{debris}+N_{agg}$ & release rates, injury state & measured particle count \\",
    r"Purity score & Ratio or classifier & $Q_{EV}=\frac{N_{bona\ fide\ EV}}{N_{meas}}$ & marker panel, debris state & EV quality score \\": r"Purity score & Ratio or classifier \cite{webber2013pure,brennan2020comparison,welsh2024minimal} & $Q_{EV}=\frac{N_{bona\ fide\ EV}}{N_{meas}}$ & marker panel, debris state & EV quality score \\",
    r"Acceptable window & Constrained optimization & maximize $Y_{EV}$ subject to $Q_{EV}>Q^*$ and $Viability>V^*$ & dose, cell state, harvest & operating region \\": r"Acceptable window & Constrained optimization \cite{orlacchio2023effects,gudvangen2022electroporation,welsh2024minimal} & maximize $Y_{EV}$ subject to $Q_{EV}>Q^*$ and $Viability>V^*$ & dose, cell state, harvest & operating region \\",
    # A15/A16
    r"Harvest kinetics & Time-dependent EV accumulation & minute-to-hour release curves after nsPEF & Determines optimal collection window & Supported \\": r"Harvest kinetics & Time-dependent EV accumulation & minute-to-hour release curves after nsPEF & Determines optimal collection window & Supported \cite{orlacchio2023effects,gudvangen2022electroporation,mendt2018generation} \\",
    r"Cell-normalized yield & EV output per cell or protein & NTA/flow/ELISA normalized to cell number & Defines productivity & Established \\": r"Cell-normalized yield & EV output per cell or protein & NTA/flow/ELISA normalized to cell number & Defines productivity & Established \cite{welsh2024minimal,paganini2019scalable,mendt2018generation} \\",
    r"Isolation recovery & SEC, ultracentrifugation, density gradient, affinity capture & recovery and purity comparisons & Converts released EVs into isolated product & Established \\": r"Isolation recovery & SEC, ultracentrifugation, density gradient, affinity capture & recovery and purity comparisons & Converts released EVs into isolated product & Established \cite{watson2018scalable,brennan2020comparison,webber2013pure} \\",
    r"QC characterization & Size, morphology, markers, contaminants & NTA, cryo-EM/TEM, western blot/flow, proteomics & Determines whether product meets EV criteria & Established \\": r"QC characterization & Size, morphology, markers, contaminants & NTA, cryo-EM/TEM, western blot/flow, proteomics & Determines whether product meets EV criteria & Established \cite{welsh2024minimal,jeppesen2019reassessment,kugeratski2021quantitative} \\",
    r"Batch variability & Donor, passage, culture, and pulse variability & repeated-batch studies, mixed-effects data & Determines reproducibility and manufacturability & Supported \\": r"Batch variability & Donor, passage, culture, and pulse variability & repeated-batch studies, mixed-effects data & Determines reproducibility and manufacturability & Supported \cite{paganini2019scalable,watson2018scalable,mendt2018generation,thakur2024global} \\",
    r"Scalable flow processing & Closed-system pulse chamber and inline monitoring & field uniformity, sterility, temperature, flow control & Enables translational manufacturing & Hypothesized \\": r"Scalable flow processing & Closed-system pulse chamber and inline monitoring & field uniformity, sterility, temperature, flow control & Enables translational manufacturing & Hypothesized \cite{liu2026microfluidic,garcia2016microfluidic,bonakdar2017microfluidic,thakur2024global} \\",
    r"Cell-normalized yield & Productivity metric & $Y_{EV}=\frac{N_{EV}}{N_{cells}}$ & EV count, cell number & EV yield per cell \\": r"Cell-normalized yield & Productivity metric \cite{welsh2024minimal,paganini2019scalable,mendt2018generation} & $Y_{EV}=\frac{N_{EV}}{N_{cells}}$ & EV count, cell number & EV yield per cell \\",
    r"Harvest accumulation & Time-integral model & $N_{EV}(t_h)=\int_0^{t_h}R_{EV}(t)dt$ & release rate, harvest time & collected EVs \\": r"Harvest accumulation & Time-integral model \cite{orlacchio2023effects,gudvangen2022electroporation,mendt2018generation} & $N_{EV}(t_h)=\int_0^{t_h}R_{EV}(t)dt$ & release rate, harvest time & collected EVs \\",
    r"Isolation recovery & Process-yield balance & $N_{iso}=\eta_{iso}N_{released}$ & isolation method, recovery & isolated EV count \\": r"Isolation recovery & Process-yield balance \cite{watson2018scalable,brennan2020comparison,webber2013pure} & $N_{iso}=\eta_{iso}N_{released}$ & isolation method, recovery & isolated EV count \\",
    r"Purity metric & Product-quality ratio & $Purity=\frac{N_{EV}}{Protein_{total}}$ or marker-based score & particle/protein/marker data & purity estimate \\": r"Purity metric & Product-quality ratio \cite{webber2013pure,brennan2020comparison,welsh2024minimal} & $Purity=\frac{N_{EV}}{Protein_{total}}$ or marker-based score & particle/protein/marker data & purity estimate \\",
    r"Batch model & Mixed-effects model & $Y_{b}=\mu+\alpha_{cell}+\beta_{dose}+u_b+\epsilon$ & cell state, dose, batch & reproducibility \\": r"Batch model & Mixed-effects model \cite{paganini2019scalable,watson2018scalable,mendt2018generation,thakur2024global} & $Y_{b}=\mu+\alpha_{cell}+\beta_{dose}+u_b+\epsilon$ & cell state, dose, batch & reproducibility \\",
    r"Optimization objective & Constrained process optimization & maximize $Potency\times Y_{iso}$ subject to purity and viability constraints & dose, harvest, isolation & optimal protocol \\": r"Optimization objective & Constrained process optimization \cite{paganini2019scalable,watson2018scalable,mendt2018generation,thakur2024global,welsh2024minimal} & maximize $Potency\times Y_{iso}$ subject to purity and viability constraints & dose, harvest, isolation & optimal protocol \\",
}


def append_missing_bib_entries() -> None:
    text = BIB.read_text(encoding="utf-8")
    additions: list[str] = []
    for block in BIB_ENTRIES.strip().split("\n\n@"):
        entry = block if block.startswith("@") else "@" + block
        key = entry.split("{", 1)[1].split(",", 1)[0]
        if f"{{{key}," not in text:
            additions.append(entry)
    if additions:
        BIB.write_text(text.rstrip() + "\n\n" + "\n\n".join(additions) + "\n", encoding="utf-8")


def patch_main() -> None:
    text = MAIN.read_text(encoding="utf-8")
    for old, new in REPLACEMENTS.items():
        count = text.count(old)
        if count != 1:
            raise RuntimeError(f"Expected one match for replacement, found {count}: {old[:100]}")
        text = text.replace(old, new)
    MAIN.write_text(text, encoding="utf-8")


def main() -> int:
    append_missing_bib_entries()
    patch_main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
