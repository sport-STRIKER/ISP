import pysam
from pyfaidx import Fasta
from collections import defaultdict

REF_NAME = "chr1"

MIN_COVERAGE = 2                
MIN_ALLELE_FREQ = 0.2

def PotentialVariants(samfile, fasta, variants, total_coverage):
    for read in samfile:
        if not read.is_unmapped and read.mapping_quality >= 20:  
            ref_seq = fasta[REF_NAME][read.reference_start:read.reference_end].seq.upper()
            query_seq = read.query_sequence

            for i, (ref_base, query_base) in enumerate(zip(ref_seq, query_seq)):
                pos = read.reference_start + i + 1  
                total_coverage[pos] += 1 

                if ref_base != query_base and ref_base != "N" and query_base != "N":
                    variants[pos][query_base] += 1   

def CheckVariants(samfile, fasta):
    variants = defaultdict(lambda: defaultdict(int))
    total_coverage = defaultdict(int)

    PotentialVariants(samfile, fasta, variants, total_coverage)

    vcf_lines = []
    for pos in variants:
        coverage = total_coverage[pos]
        if coverage >= MIN_COVERAGE:
            for alt_base, count in variants[pos].items():
                freq = count / coverage
                if count >= MIN_COVERAGE and freq > MIN_ALLELE_FREQ:  
                    ref_base = fasta[REF_NAME][pos-1:pos].seq.upper() 
                    vcf_lines.append(f"{REF_NAME}\t{pos}\t.\t{ref_base}\t{alt_base}")

    return vcf_lines

def PrintResult(vcf_lines):
    with open("output.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\t\t\tID\tREF\tALT\n")
        
        for line in sorted(vcf_lines):
            print(line, file=f)

    print(f"Варианты сохранены в output.vcf ({len(vcf_lines)} SNP)")

def main():
    fasta = Fasta("chr1.fasta")
    samfile = pysam.AlignmentFile("test_no_indel.sam", "r")

    PrintResult(CheckVariants(samfile, fasta))

    samfile.close()
    fasta.close()

main()
