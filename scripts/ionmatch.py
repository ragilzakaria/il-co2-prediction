"""Identify ions in xyz geometries by matching heavy-atom+H connectivity graphs to our ion list."""
import csv
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
RDLogger.DisableLog('rdApp.*')
IONS = str(Path(__file__).resolve().parents[1] / 'data' / 'ions' / 'ions_with_smiles.csv')

def skeleton(m):
    m = Chem.RWMol(m)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    for a in m.GetAtoms():
        a.SetFormalCharge(0); a.SetIsAromatic(False); a.SetNoImplicit(True); a.SetNumExplicitHs(0); a.SetNumRadicalElectrons(0)
    m = m.GetMol(); m.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(m, canonical=True, allHsExplicit=False)

def ion_table():
    t = {}
    for r in csv.DictReader(open(IONS)):
        m = Chem.MolFromSmiles(r['smiles'])
        if m is None: continue
        t[skeleton(Chem.AddHs(m))] = r['ion_abbreviation']
    return t

def frames(path):
    lines = open(path).read().splitlines(); i = 0
    while i < len(lines):
        if not lines[i].strip(): i += 1; continue
        n = int(lines[i].split()[0]); title = lines[i+1]
        at=[]
        for l in lines[i+2:i+2+n]:
            el,x,y,z=l.split()[:4]
            el=el.strip().capitalize()
            at.append(f"{el} {float(x):.6f} {float(y):.6f} {float(z):.6f}")
        yield title, '\n'.join([str(n), ''] + at)
        i += 2 + n

def fragments(xyzblock):
    m = Chem.MolFromXYZBlock(xyzblock)
    rdDetermineBonds.DetermineConnectivity(m)
    out = []
    for f in Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False):
        f.UpdatePropertyCache(strict=False)
        out.append((skeleton(f), CalcMolFormula(f), f.GetNumAtoms()))
    return out
