class InternalRBSChecker:
    """
    Detects potential internal ribosome binding sites (RBS) inside a CDS.
    A simple heuristic looks for a Shine-Dalgarno-like motif followed by
    a start codon within a short spacer distance.
    """

    def __init__(self):
        self.sd_motifs = ["AGGAGG", "GGAGG", "AGGA", "GAGG"]
        self.start_codons = ["ATG", "GTG", "TTG"]
        self.min_spacer = 4
        self.max_spacer = 12

    def initiate(self) -> None:
        pass

    def run(self, seq: str):
        """
        Returns:
            (True, None) if no internal RBS is found
            (False, problematic_region) if one is found
        """
        seq = seq.upper()

        for motif in self.sd_motifs:
            start = 0
            while True:
                i = seq.find(motif, start)
                if i == -1:
                    break

                motif_end = i + len(motif)

                for spacer in range(self.min_spacer, self.max_spacer + 1):
                    codon_start = motif_end + spacer
                    codon_end = codon_start + 3

                    if codon_end <= len(seq):
                        codon = seq[codon_start:codon_end]
                        if codon in self.start_codons:
                            bad_region = seq[i:codon_end]
                            return False, bad_region

                start = i + 1

        return True, None