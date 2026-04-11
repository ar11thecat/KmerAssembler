from pathlib import Path
from random import randrange

import alg_utils as alg
from sequence import Sequence


def main():

    file_path = Path(__file__).parent.parent / "data" / "aromatase_aa.txt"
    seq_string = file_path.read_text().replace("\n", "")

    k = 6
    
    my_seq = Sequence(seq_string, "protein")
    print(f"Starting with sequence:\n{my_seq}\n")
    print(f"Length: {len(my_seq)}\n")

    my_kmer_dict = my_seq.kmer_dict(k)
    print(f"kmer dictionary:\n{my_kmer_dict}\n")

    my_graph = alg.kmer_dict2graph(my_kmer_dict)
    print(f"kmer graph (nodes = k-1mers, edges = kmers):\n{my_graph}\n")

    my_path = alg.eulerian_path(my_graph)
    print(f"Eulerian path:\n{my_path}\n")

    my_assembled = Sequence(alg.condense(my_path), "protein")
    print(f"Assembled sequence:\n{my_assembled}\n")
    print(f"The assembled sequence is the same as the original: {my_assembled == my_seq}\n")

    return

    all_paths = alg.all_eulerian_paths(my_graph)

    print(f"Possible assembled sequences:")
    for path in all_paths:
        assembled_path = Sequence(alg.condense(path), "protein")
        print(assembled_path)

    return

if __name__ == "__main__":
    print("Hello World!\n")
    main()
