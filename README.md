```
                   /$$                  /$$$$$$$                    /$$     /$$$$$$          
                  | $$                 | $$__  $$                  | $$    /$$__  $$         
  /$$$$$$ /$$   /$| $$ /$$$$$$  /$$$$$$| $$  \ $$/$$$$$$  /$$$$$$ /$$$$$$ | $$  \ $$ /$$$$$$$
 /$$__  $| $$  | $| $$/$$__  $$/$$__  $| $$$$$$$/$$__  $$/$$__  $|_  $$_/ | $$$$$$$$/$$_____/
| $$  \ $| $$  | $| $| $$$$$$$| $$  \__| $$____| $$$$$$$| $$  \ $$ | $$   | $$__  $|  $$$$$$ 
| $$  | $| $$  | $| $| $$_____| $$     | $$    | $$_____| $$  | $$ | $$ /$| $$  | $$\____  $$
|  $$$$$$|  $$$$$$| $|  $$$$$$| $$     | $$    |  $$$$$$| $$$$$$$/ |  $$$$| $$  | $$/$$$$$$$/
 \______/ \____  $|__/\_______|__/     |__/     \_______| $$____/   \___/ |__/  |__|_______/ 
          /$$  | $$                                     | $$                                 
         |  $$$$$$/                                     | $$                                 
          \______/                                      |__/                                           
```
# Welcome to oylerPeptideAssembler!
> [!WARNING]
> ⌛ sorry! Work in progress ⌛

## PEAKS Simulator
The PEAKS Simulator simulates a Peptide-Spectrum Matching (PSM) report with the format produced by the PEAKS software.
It is used to simulate realistic peptides with `m/z`, `z`, `mass`, and `local confidence` values, to be used in the downstream assembly.

### Input / Output
Takes as input `.txt` file with the aminoacid sequence of the target protein, and a `.txt` file with the aminoacid sequence to use as noise.\
It can produce the results as `.csv` file or as `pandas DataFrame`.

### What it can do:
- [x] Simulate the cuts by a selectable range of proteases (tripsin, chymotripsin, elastase and others), considering cut probabilities and proline suppression
- [x] Consider the protonation of the residues to derive the m/z and z measures, as it would happen in MS1
- [x] Correct the neutral mass of the peptide by removing tryptic water
- [x] Simulate the local confidence by considering:
+ the intrinsic confidence of each residue (for example, low for Leucine and Isoleucine, which have the same mass)
+ the position of the residue (residues close to the N-terminus of C-terminus can only form shorter b and y ions, which are harder to detect)
+ gaussian noise

### What it cannot do:
- [ ] PTMs are not yet modelled\
- [ ] Only the 20 canonycal aminoacid letters are allowed\
- [ ] The protonation count is approximated, and doesn't model precise chemo-physical interactions\
- [ ] The local confidence model is an approximation, and is not derived from spectral reconstruction\
