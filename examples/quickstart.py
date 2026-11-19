"""Detect pockets on 1KNE (chromodomain + methylated histone peptide) and score them
against the bound peptide.   python examples/quickstart.py"""
from pathlib import Path

import ppdpy

pdb = Path(__file__).with_name("1KNE.pdb")
res = ppdpy.detect_pockets(str(pdb), verbose=True)   # chain A = receptor, short chain P = peptide (removed)

print(res)
for row in res.evaluate():
    print(f"pocket {row['pocket']}: residue recall {row['residue_recall']:.2f}, "
          f"precision {row['residue_precision']:.2f}, DCC {row['dcc']:.1f} A")

best = res.by_volume()[0]
print("lining residues of the largest pocket:", " ".join(best.residues))
res.save("1KNE_pockets")     # json / csv / pdb (pocket id in B-factor) / PyMOL script / log
