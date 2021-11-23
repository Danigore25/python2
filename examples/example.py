from Bio.Seq import Seq
from Bio.SeqUtils import GC

print("Write the DNA sequence: ")
dna = input()

print("5'", Seq.complement, "3'")
print("3'", Seq.reverse_complement, "5'")

print('GC Content: ', GC(dna))

rna = dna.translate
