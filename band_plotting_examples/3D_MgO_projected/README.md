# MgO — band structure and projected (fat) bands

Sample output of `$pdosper bands on` for the **Band Structure Plotter** page of RIPER-Tools:
upload `bands.xyz` (plain bands) or `bands_proj.xyz` (fat bands) there. Both files come from
the same run (`output`), whose `control` is included.

MgO rocksalt, conventional cubic cell (4 Mg + 4 O), PBE / pob-DZVP-rev2, 4x4x4 k-mesh,
path Γ-X-M-Γ (61 points). The relevant `control` groups:

```
$kpoints
  nkpoints 4 4 4
  kptlines 3
  recipr 0.0 0.0 0.0  0.5 0.0 0.0  20
  recipr 0.5 0.0 0.0  0.5 0.5 0.0  20
  recipr 0.5 0.5 0.0  0.0 0.0 0.0  21
$pdosper
  elements mg,o
  bands on
```

`bands_proj.xyz` therefore carries 8 weight buckets (Mg/O × total,s,p,d), listed in its
header. Useful numbers for the plotter: Fermi level (gap middle) = 0.3840051400 a.u.;
the VBM is band 40 at 0.291331 a.u. (the automatic detection of the plotter finds this,
gap 5.044 eV, direct at Γ). The valence band is ~99% O-2p, the conduction band ~92% Mg
(s at Γ, p towards M) — watch the lowest conduction band change color in the fat-band
plot. The command-line alternative to the web page is included:

```bash
python3 bands_plotter_adapted.py normal    bands.xyz      0.3840051400 61 -8 12 "MgO bands (PBE)"     mgo_bands.png "Γ,X,M,Γ"
python3 bands_plotter_adapted.py projected bands_proj.xyz 0.3840051400 61 -8 12 "MgO fat bands (PBE)" mgo_proj.png  "Γ,X,M,Γ" "element o  total" "element mg  total"
```

Note that the weights of individual bands inside a degenerate multiplet depend on the
arbitrary eigenvector basis within the multiplet; only their sum over the multiplet is well
defined. Fat-band plots show the multiplet together, so they are unaffected.
