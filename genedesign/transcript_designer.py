import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.aminoAcidToCodons = {}
        self.stopCodons = ["TAA", "TGA", "TAG"]

        self.rbsChooser = None

        self.forbiddenChecker = None
        self.hairpinChecker = None
        self.codonChecker = None
        self.promoterChecker = None
        self.internalRBSChecker = None

        self.random_seed = 134
        self.num_candidates = 75

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        random.seed(self.random_seed)

        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # All synonymous codons for E. coli design
        self.aminoAcidToCodons = {
            'A': ["GCT", "GCC", "GCA", "GCG"],
            'C': ["TGT", "TGC"],
            'D': ["GAT", "GAC"],
            'E': ["GAA", "GAG"],
            'F': ["TTT", "TTC"],
            'G': ["GGT", "GGC", "GGA", "GGG"],
            'H': ["CAT", "CAC"],
            'I': ["ATT", "ATC", "ATA"],
            'K': ["AAA", "AAG"],
            'L': ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
            'M': ["ATG"],
            'N': ["AAT", "AAC"],
            'P': ["CCT", "CCC", "CCA", "CCG"],
            'Q': ["CAA", "CAG"],
            'R': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            'S': ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
            'T': ["ACT", "ACC", "ACA", "ACG"],
            'V': ["GTT", "GTC", "GTA", "GTG"],
            'W': ["TGG"],
            'Y': ["TAT", "TAC"]
        }

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.codonChecker = CodonChecker()
        self.promoterChecker = PromoterChecker()
        self.internalRBSChecker = InternalRBSChecker()

        if hasattr(self.forbiddenChecker, "initiate"):
            self.forbiddenChecker.initiate()
        if hasattr(self.hairpinChecker, "initiate"):
            self.hairpinChecker.initiate()
        if hasattr(self.codonChecker, "initiate"):
            self.codonChecker.initiate()
        if hasattr(self.promoterChecker, "initiate"):
            self.promoterChecker.initiate()
        if hasattr(self.internalRBSChecker, "initiate"):
            self.internalRBSChecker.initiate()

    def _generate_candidate_codons(self, peptide: str) -> list[str]:
        #generates on synonymous codon sequence for the peptide and appends a stop codon

        codons = []
        for aa in peptide:
            if aa not in self.aminoAcidToCodons:
                raise ValueError(f"Unsupported amino acid: {aa}")
            
            choices = self.aminoAcidToCodons[aa]
            weights = [self.codonChecker.codon_frequencies.get(codon, 0.01) for codon in choices]
            chosen = random.choices(choices, weights=weights, k=1)[0]
            codons.append(chosen)

        codons.append("TAA")
        return codons

    def _score_candidate(self, codons: list[str], cds: str) -> float:
        score = 0.0

        forbidden_ok, _ = self.forbiddenChecker.run(cds)
        if not forbidden_ok:
            return float("-inf")

        promoter_ok, _ = self.promoterChecker.run(cds)
        if promoter_ok:
            score += 2.0
        else:
            score -= 8.0

        internal_rbs_ok, _ = self.internalRBSChecker.run(cds)
        if internal_rbs_ok:
            score += 2.0
        else:
            score -= 8.0

        _, hairpins = hairpin_checker(cds)
        hairpin_count = 0 if hairpins is None else len(hairpins)
        score -= 1.0 * hairpin_count

        codon_ok, codon_diversity, rare_codon_count, cai_value = self.codonChecker.run(codons)

        if codon_ok:
            score += 3.0
        else:
            score -= 3.0

        score += 10.0 * cai_value
        score += 1.0 * codon_diversity
        score -= 1.0 * rare_codon_count

        return score

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Generates multiple CDS candidates, evaluates them, and returns
        the best transcript. Prefers fully valid candidates, but falls
        back to the best available candidate instead of raising an error.
        """
        best_valid_codons = None
        best_valid_cds = None
        best_valid_score = float("-inf")

        best_any_codons = None
        best_any_cds = None
        best_any_score = float("-inf")

        for _ in range(self.num_candidates):
            codons = self._generate_candidate_codons(peptide)
            cds = ''.join(codons)
            score = self._score_candidate(codons, cds)

            if best_any_codons is None or score > best_any_score:
                best_any_score = score
                best_any_codons = codons
                best_any_cds = cds

            if score != float("-inf") and score > best_valid_score:
                best_valid_score = score
                best_valid_codons = codons
                best_valid_cds = cds

        if best_valid_codons is not None:
            selectedRBS = self.rbsChooser.run(best_valid_cds, ignores)
            return Transcript(selectedRBS, peptide, best_valid_codons)

        if best_any_codons is not None:
            selectedRBS = self.rbsChooser.run(best_any_cds, ignores)
            return Transcript(selectedRBS, peptide, best_any_codons)

        raise Exception("Unable to generate any CDS candidate.")
