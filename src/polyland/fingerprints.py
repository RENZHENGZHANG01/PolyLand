"""RDKit fingerprint helpers shared by predictive and ICL workflows."""

from __future__ import annotations

import numpy as np


def smiles_to_fingerprints(
    smiles: list[str], method: str = "Morgan", radius: int = 2, n_bits: int = 2048
) -> np.ndarray:
    """Convert polymer repeat-unit SMILES to a dense fingerprint matrix."""

    try:
        from rdkit import Chem
        from rdkit.Chem import MACCSkeys, rdMolDescriptors
    except ImportError as exc:  # pragma: no cover - optional scientific dependency
        raise RuntimeError("Install PolyLand with the 'ml' extra to use fingerprints") from exc

    values: list[np.ndarray] = []
    for value in smiles:
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {value}")
        if method == "Morgan":
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        elif method == "RDKit":
            fp = Chem.RDKFingerprint(mol, maxPath=radius, fpSize=n_bits)
        elif method == "MACCS":
            fp = MACCSkeys.GenMACCSKeys(mol)
        elif method == "TopologicalTorsion":
            fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
                mol, nBits=n_bits
            )
        elif method == "AtomPair":
            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
        else:
            raise ValueError(f"Unsupported fingerprint method: {method}")
        values.append(np.asarray(fp, dtype=np.float32))
    return np.vstack(values)

