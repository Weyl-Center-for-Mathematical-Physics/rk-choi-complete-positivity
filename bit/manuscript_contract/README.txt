Journal of Computational Physics source package, version 3.5.

Primary files:
- Manuscript_JCP.tex: line-numbered Elsevier review manuscript for submission
- Manuscript_JCP_clean.tex: two-column author convenience build
- Supplementary_Material_JCP.tex: derivations and exact certificates
- Cover_Letter_JCP.tex: journal-specific cover letter
- Highlights.txt: separate submission highlights
- Figure_Captions.txt: captions synchronized with the article source
- references_jcp.tex: authoritative embedded bibliography
- references.bib: machine-readable companion bibliography for added root/software citations
- figures/: vector article figures plus the graphical abstract

Version 3.5 is an editorial-polish release built on the unchanged v3.4 mathematical and computational results. It tightens the abstract, introduction, step-control discussion, conclusion, Supplementary Material, cover letter, captions, and graphical abstract while preserving the theorem scopes and fail-closed certification rules. The work proxy counts RK stages and exponential actions; algebraic certification costs are reported separately.

The graphical abstract is generated deterministically by Matplotlib from certified data and equations. Author order, affiliations, addresses, emails, ORCIDs, CRediT roles, the competing-interest statement, originality, author approval, and the absence of concurrent review have been confirmed. The Zenodo DOI placeholder must be replaced before submission.

Reproducibility notes:
- PDF is the authoritative vector figure format; EPS files flatten partial transparency.
- The staged reproduction workflow and pinned environment are documented in the code-and-data archive.
