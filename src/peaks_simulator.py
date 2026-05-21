from collections import Counter
import math
from pandas import DataFrame
from pathlib import Path
import random
import re
from typing import Callable


PROTON_MASS = 1.00728

WATER_MASS = 18.01056

AA_MASS = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}


class PeaksSimulator:

    attributes = (
        "Scan", "Peptide", "Tag Length", "ALC (%)", "length", "m/z", "z", "RT",
        "Predict RT", "Area", "Mass", "ppm", "PTM", "local confidence (%)", "tag (>=0%)", "mode",
    )

    residue_bias = {
        'A': 90, 'R': 95, 'N': 85, 'D': 90,
        'C': 90, 'E': 90, 'Q': 85, 'G': 95,
        'H': 90, 'I': 75, 'L': 75, 'K': 95,
        'M': 85, 'F': 90, 'P': 95, 'S': 90,
        'T': 90, 'W': 90, 'Y': 90, 'V': 90,
    }

    proteases = {
        "tripsin": {"primary_sites": "KR", "primary_prob": 0.95, "proline_prob": 0.05},
        "chymotrypsin": {"primary_sites": "FYW", "primary_prob": 0.80, "proline_prob": 0.00},
        "lys_c": {"primary_sites": "K", "primary_prob": 0.98},
        "glu_c": {"primary_sites": "E", "primary_prob": 0.85},
        "elastase": {"primary_sites": "AVGSLI", "primary_prob": 0.35},
    }


    @staticmethod
    def _protease_factory(primary_sites: str,
                          primary_prob: float,
                          proline_prob: float=None) -> Callable:
        # Generates a function that simulates the cuts of a peptidase given parameters
        
        query = f"([{primary_sites}])(P?)"

        def protease(seq: str) -> tuple[list[int], list[int]]:
            cut_indices = []
            cut_prob = []

            for match in re.finditer(query, seq):
                cut_index = match.end(1)
                has_proline = bool(match.group(2))
                use_proline = has_proline and (proline_prob is not None)
                chance_to_cut = proline_prob if use_proline else primary_prob

                cut_indices.append(cut_index)
                cut_prob.append(chance_to_cut)

            return cut_indices, cut_prob

        return protease


    def __init__(self, dir: str):
        self.directory = Path(__file__).parent.parent / dir


    def simulate_df(self,
                      target_prot: str,
                      noise_prot: str,
                      fold: int,
                      enrichment: float,
                      proteases: list[str],
                      scans: int) -> DataFrame:
        simulated = self._simulate(target_prot, noise_prot, fold, enrichment, proteases, scans)
        return DataFrame.from_dict(simulated)


    def simulate_file(self,
                       target_prot: str,
                       noise_prot: str,
                       fold: int,
                       enrichment: float,
                       proteases: list[str],
                       scans: int):
        simulated_file = self.directory / f"{target_prot}.sim-peaks.csv"
        
        simulated = self._simulate(target_prot, noise_prot, fold, enrichment, proteases, scans)
        simulated_df = DataFrame.from_dict(simulated)
        simulated_df.to_csv(simulated_file, header=True, index=False)
        

    def _simulate(self,
                  target_prot: str,
                  noise_prot: str,
                  fold: int, enrichment: float,
                  proteases: list[str],
                  scans: int) -> dict:
        
        prot_file = self.directory / f"{target_prot}.txt"
        noise_file = self.directory / f"{noise_prot}.txt"
        prot_seq = prot_file.read_text().strip()
        noise_seq = noise_file.read_text().strip()

        peptides = self._cut(prot_seq, proteases, math.ceil(fold * enrichment))
        noise = self._cut(noise_seq, proteases, math.ceil(fold * (1-enrichment)))
        peptides += noise
                
        sampled_peptides = self._sample(peptides, scans)

        simulated = {attribute: [] for attribute in PeaksSimulator.attributes}
        
        for i, peptide in enumerate(sampled_peptides):
            m_z, z, mass = self._acquire(peptide)
            local_confidence = self._local_confidence(peptide)

            simulated["Scan"].append(f"F1:{i:04d}")
            simulated["Peptide"].append(peptide)
            simulated["Tag Length"].append(len(peptide))
            simulated["ALC (%)"].append(round(sum(local_confidence) / len(peptide), 4))
            simulated["length"].append(len(peptide))
            simulated["m/z"].append(round(m_z, 4))
            simulated["z"].append(z)
            simulated["RT"].append('')
            simulated["Predict RT"].append('')
            simulated["Area"].append('')
            simulated["Mass"].append(round(mass, 4))
            simulated["ppm"].append('')
            simulated["PTM"].append('')
            simulated["local confidence (%)"].append(' '.join(map(str, local_confidence)))
            simulated["tag (>=0%)"].append('')
            simulated["mode"].append('')

        return simulated
        

    def _cut(self, seq: str, proteases: list[str], fold: int) -> list[str]:
        # cut the sequence simulating the combined effect of the proteases
        
        peptides = []
        cut_indices = []
        cut_probs = []
        
        for p in proteases:
            protease = PeaksSimulator._protease_factory(**PeaksSimulator.proteases[p])
            indices, probs = protease(seq)
            cut_indices.extend(indices)
            cut_probs.extend(probs)

        for _ in range(fold):
            indices = []
            for i, p in zip(cut_indices, cut_probs):
                if random.random() < p:
                    indices.append(i)
            indices = [0] + sorted(indices) + [len(seq)]
            peptides.extend([seq[indices[i-1]:indices[i]] for i in range(1, len(indices))
                            if indices[i] > indices[i-1]])

        return peptides


    def _sample(self, raw_peptides: list[str], scans: int) -> list[str]:
        # simulate the sampling of LC columns and Mass Spectroscopy

        filtered_peptides = [p for p in raw_peptides if 6 <= len(p) <= 35]

        peptide_counts = Counter(filtered_peptides)
        detectability = []

        for peptide, count in peptide_counts.items():
            basic_residues = len(re.findall("([RKH])", peptide))
            ionizability_score = 1.5 if basic_residues > 0 else 0.2
            detectability.append(count * ionizability_score)

        sampled_peptides = random.choices(
            population=list(peptide_counts.keys()),
            weights=detectability,
            k=scans,
        )

        return sampled_peptides


    def _acquire(self, peptide: str) -> tuple[float, int, float]:
        # simulate the m/z and z acquisition and mass calculation
        
        basic_residues = len(re.findall("([RKH])", peptide)) + 1 # includes N-terminal
        z = random.randint(1, min(basic_residues, 3)) # assumes z > 0 as _sample already filtered by ionizability

        mass = sum([AA_MASS[aa] for aa in peptide]) + WATER_MASS
        m_z = (mass + z * PROTON_MASS) / z

        return m_z, z, mass


    def _local_confidence(self, peptide: str) -> list[int]:
        # simulate the local confidence evaluation
        
        local_confidence = []

        for i, aa in enumerate(peptide):
            residue_bias = PeaksSimulator.residue_bias[aa]
            position_penalty = -int(10 * math.exp(-0.8 * i) + 5 * math.exp(-0.8 * (len(peptide)-i-1)))
            noise = int(random.gauss(0, 3))
            
            local_confidence.append(residue_bias + position_penalty + noise)

        return local_confidence
