# LiF — band structure and projected (fat) bands

Sample output of `$pdosper bands on` for the **Band Structure Plotter** page of RIPER-Tools:
upload `bands.xyz` (plain bands) or `bands_proj.xyz` (fat bands) there. Both files come from
the same run (`output`), whose `control` is included.

LiF rocksalt, conventional cubic cell (4 Li + 4 F), PBE / pob-DZVP-rev2, 4x4x4 k-mesh,
path Γ-X-M-Γ (61 points). Here `$pdosper` was used **without** any selection:

```
$pdosper
  bands on
```

so `bands_proj.xyz` carries the full default set of 40 weight buckets — every atom and every
element × (total, s, p, d) — a good file to explore the bucket picker of the plotter page.
Useful numbers: Fermi level (gap middle) = 0.2009489415 a.u.; the VBM is band 24 at
0.031130 a.u. (automatic detection: gap 9.242 eV, direct at Γ). The valence band is ~99% F-2p
with ~5% Li admixture; the conduction band minimum at Γ is 97% Li. Command line:

```bash
python3 bands_plotter_adapted.py normal    bands.xyz      0.2009489415 61 -10 16 "LiF bands (PBE)"     lif_bands.png "Γ,X,M,Γ"
python3 bands_plotter_adapted.py projected bands_proj.xyz 0.2009489415 61 -10 16 "LiF fat bands (PBE)" lif_proj.png  "Γ,X,M,Γ" "element f  total" "element li  total"
```

Note that the weights of individual bands inside a degenerate multiplet depend on the
arbitrary eigenvector basis within the multiplet; only their sum over the multiplet is well
defined. Fat-band plots show the multiplet together, so they are unaffected.
