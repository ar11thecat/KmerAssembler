import decors as dec


class Sequence:

    def __init__(self, sequence: str, category: str):
        self._alphabet = None
        self._sequence = None
        self.alphabet = category
        self.sequence = sequence


    @property
    def alphabet(self) -> str:
        return self._alphabet
    

    @property
    def sequence(self) -> str:
        return self._sequence


    @alphabet.setter
    def alphabet(self, category: str) -> None:
        match category:
            case "dna": alph = "ATCG"
            case "protein": alph = "ACDEFGHIKLMNPQRSTVWY"
            case _: raise ValueError("Unsupported sequence type")
            
        self._alphabet = set(alph)
    

    @sequence.setter
    def sequence(self, sequence: str) -> None:
        if not sequence:
            raise ValueError("Sequence cannot be empty")

        norm_sequence = sequence.upper().strip()

        if not all(base in set(self._alphabet) for base in norm_sequence):
            invalid = set(norm_sequence) - set(self._alphabet)
            raise ValueError(f"Invalid letters in sequence: {invalid}")

        self._sequence = norm_sequence


    def __str__(self) -> str:
        return self._sequence


    def __len__(self) -> int:
        return len(self._sequence)


    def __eq__(self, other) -> bool:
        if isinstance(other, Sequence):
            return self._sequence == other._sequence

        return False


    def __add__(self, other):
        if isinstance(other, Sequence):
            return Sequence(self._sequence + other._sequence)

        raise TypeError(f"Cannot add Sequence with {type(other)}")


    @dec.timer
    def _all_kmers(self, k: int) -> list:
        kmers = list(self._alphabet)

        def add_letter(imers: list) -> list:
            jmers = []
            for imer in imers:
                for letter in self._alphabet:
                    jmers.append(imer + letter)
            return jmers

        for i in range(k - 1):
            kmers = add_letter(kmers)

        return kmers


    @dec.timer
    def all_kmer_dict(self, k: int) -> dict:
        all_kmers = self._all_kmers(k)

        kmers = dict.fromkeys(all_kmers, 0)

        for i in range(len(self) - k + 1):
            kmer = self._sequence[i:i+k]
            kmers[kmer] += 1

        return kmers


    @dec.timer
    def kmer_dict(self, k: int) -> dict:
        kmers = dict()

        for i in range(len(self) - k + 1):
            kmer = self._sequence[i:i+k]
            kmers.setdefault(kmer, 0)
            kmers[kmer] += 1

        return kmers
