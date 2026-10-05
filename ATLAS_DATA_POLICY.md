# Atlas data policy

The public release candidate contains no atlas NIfTI images and no derived atlas operator matrices.

For reproducibility it may publish:
- frozen SHA-256 identities;
- matrix dimensions and nonzero counts;
- construction/resampling recipes;
- source citations and upstream license notes;
- user-local validation tools.

It must not redistribute the five frozen project atlas image bytes unless exact upstream provenance and redistribution rights are established for those exact bytes/derivatives.

Current conservative disposition:
- AAL116 legacy project byte: no redistribution;
- Dosenbach161 derived project byte: no redistribution;
- EZ116/Juelich-family legacy project byte: no redistribution;
- Harvard-Oxford111 project byte: no redistribution in this release despite current FSL Harvard-Oxford CC BY-SA 4.0, because the exact project derivative chain is not yet frozen;
- Talairach-Tournoux97 project byte: no redistribution.
