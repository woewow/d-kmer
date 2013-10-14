import numpy as np
import os, sys, re
from cPickle import load, dump
import itertools
import collections
import subprocess
import pathsetup

t200mutationlist_file = "/workspace/zchong/work/Suspicious/T200MutationOrVariantSomaticReport.txt"

standard_codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', }
stop_codons  = [ 'TAA', 'TAG', 'TGA', ]
start_codons = [ 'TTG', 'CTG', 'ATG', ]


def registerset(d, key, value):
    if key in d:
        d[key].add(value)
    else:
        d[key] = set([value])

def registerlist(d, key, value):
    if key in d:
        d[key].append(value)
    else:
        d[key] = [value]

kmer_len = 51
class Mutation():
    
    def __init__(self, fields):

        self.sample_ID, self.Allele_Freq, self.coverage, self.HGNC, self.HGNC_AAS, self.Uploaded_variation, self.CDS_position, self.Protein_position, self.Amino_acids, self.db_SNP, self.Codons, self.Consequence, self.Condel, self.Poly_Phen, self.SIFT, self.Mutation_Assessor, self.QUAL, self.COSMIC_EXACT_MATCH, self.COSMIC_LOCATION_EXACT_MATCH, self.COSMIC_LOCATION_COVERED_MATCH, self.COSMIC_LOCATION_RANGE_COVERED_MATCH, self.Allele, self.Gene, self.Feature, self.Feature_type, self.c_DNA_position, self.Existing_variation, self.HGVSc, self.HGVSp, self.DOMAINS, self.CCDS, self.ENSP, self.Ref_Seq, self.INTRON, self.EXON, self.CANONICAL, self.id, self.strand, self.Paired, self.Mutation_Type, self.MasterLog_ID, self.Mapped_SAccession, self.M_Accession, self.Mapped_Lib_ID, self.MRN, self.Imported_Date, self.Create_By, self.Update_Date, self.Update_By, self.Cut_Off_Freq, self.Can_Dr_A_CTS, self.Can_Dr_A_GEN, self.Result_Matches_MDL, self.Sample_Matches_MDL, self.Reviewed, self.Skip_Reporting, self.Comments, self.genomes1000_allele_frequency, self.Flowcell, self.Clinic, self.Disease_Type = fields

        self.chrm, self.loc, self.ref, self.var = re.split('_|/', self.Uploaded_variation)
        self.loc  = int(self.loc)
        self.Protein_position = int(self.Protein_position)
        self.chrm = "chr"+self.chrm

    def __repr__(self):

        return "<Mutation: %s>" % (self.HGNC_AAS, )

def parse_fasta(fa_filename):

    with open(fa_filename) as f:
        for line in f:
            if line.startswith("<"):
                continue
            else:
                yield line.strip()

# patient data
# # chrm2hotspots = {}
# prot2muts = {}
# with open(os.path.join(pathsetup.sharedir, 'T200MutationOrVariantSomaticReport.txt')) as f:
#     f.readline()
#     for line in f:
#         mut = Mutation(line.strip().split('\t'))
#         registerlist(prot2muts, (mut.HGNC, mut.Protein_position), mut)
#         # registerset(chrm2hotspots, mut.chrm, mut.loc)
class Contig:

    def __init__(self, chrm, start, end):

        self.chrm = chrm
        self.start = start
        self.end = end

    def __repr__(self):
        
        return "<Contig %s:%d-%d>" % (self.chrm, self.start, self.end)

def complement(base):

    return {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
    }[base]


def reverse_complement(seq):
    
    return ''.join([complement(base) for base in reversed(seq)])

class Codon(Contig):
    
    alphabet = ['A', 'T', 'G', 'C']
    
    def __init__(self, chrm, start, end, strand):

        Contig.__init__(self, chrm, start, end)
        
        self.strand = strand
        assert start+2 == end

    def __repr__(self):

        return '<Codon %s:%d(%s)>' % (self.chrm, self.start, self.strand)

    def __eq__(self, other):

        return ((self.chrm, self.strand, self.start, self.end) == (other.chrm, other.strand, other.start, other.end))

    def __hash__(self):

        return hash((self.chrm, self.start))

    def natural_seq(self):
        
        if self.strand == '-':
            return reverse_complement(self.seq)
        else:
            return self.seq

    def tr_aa(self):
        
        return standard_codon_table[self.natural_seq()]

    def snv_gen(self, rel_pos, exclude_synonymous=False, exclude_nonsense=False):
        """ 
        iterate through the single nucleotide variants
        and report whether the mutation is synonymous or not

        rel_pos: relative position (0-based)
        on the codon (0,1,2) to mutate
        
        """

        ref_aa = self.tr_aa()
        for var in Codon.alphabet: # all ATGC
            if var == self.seq[rel_pos]: # exclude reference
                continue
                
            mut = SingleMutator(self.start + rel_pos, self.seq[rel_pos], var)
            mut.ref_aa = ref_aa
            mut_codon_seql = list(self.seq)
            mut_codon_seql[rel_pos] = var
            mut_codon_seq = ''.join(mut_codon_seql)

            if self.strand == '+':
                natural_mut_seq = mut_codon_seq
            else:
                natural_mut_seq = reverse_complement(mut_codon_seq)
            if natural_mut_seq in standard_codon_table:
                mut.var_aa = standard_codon_table[natural_mut_seq]
            else:
                mut.var_aa = None

            # mut.is_synonymous = (mut.var_aa == mut.ref_aa)
            if exclude_synonymous and mut.var_aa == mut.ref_aa:
                continue

            if exclude_nonsense and not mut.var_aa:
                continue

            yield mut


class SingleMutator():

    def __init__(self, loc, ref, var):
        """ other attr:
        self.ref_aa, self.var_aa
        self.is_synonymous
        """
        
        self.loc = loc
        self.ref = ref
        self.var = var

    def __repr__(self):
        return "<Single Mutator: %d, %s->%s>" % (self.loc, self.ref, self.var)

    def mutate(self, seq, start):
        """ the seq to be mutated starts from start """

        seql = list(seq)
        assert seql[self.loc - start] == self.ref
        seql[self.loc - start] = self.var
        
        return ''.join(seql)
        
class Kmer():

    def __init__(self, seq, mutations=[]):
        """ init a kmer from sequence 
        """
        self.seq = seq
        self.mutations = mutations

def combinations(items):
    return ( set(itertools.compress(items,mask)) for mask in itertools.product(*[[0,1]]*len(items)) )

class CodonContig(Contig):

    def __init__(self, codon):
        """ initialize from one codon
        """
        Contig.__init__(self, codon.chrm, codon.start, codon.end)
        self.codons = [codon]

    def append(self, codon):
        """ append codon that is adjacent and on the right
        """
        self.end = max(self.end, codon.end)
        self.codons.append(codon)

    def kmer_gen(self, kmer_len, max_mut=2):
        """ generate kmer from the codon contig mutating
        every kmer
        """

        # shift the window from left to right
        for first in xrange(self.start, self.end - kmer_len + 2):
            last = first + kmer_len - 1

            codons_covered = [codon for codon in self.codons if codon.start <= last and codon.end >= first]
            for num_mut in xrange(1, max_mut+1):
                for codons_to_mut in itertools.combinations(codons_covered, num_mut):
                    codon_muts_list = []
                    for codon in codons_to_mut:
                        codon_muts = []     # the mutator of this codon
                        for i in xrange(3): # three bases
                            if codon.start + i >= first and codon.start + i <= last: # make sure that mutated position is within the kmer
                                codon_muts.extend(list(codon.snv_gen(i, exclude_synonymous=True, exclude_nonsense=False)))
                            
                        codon_muts_list.append(codon_muts)

                    for muts in itertools.product(*codon_muts_list):
                        kmer = self.seq[first-self.start:last+1-self.start]
                        for mut in muts:
                            kmer = mut.mutate(kmer, first)

                        yield (kmer, muts)

        # # shift the window from left to right
        # for first in xrange(self.start, self.end - kmer_len + 2):
        #     last = first + kmer_len - 1

        #     # collect all the possible mutations that can
        #     # occur to each codon
        #     codon_muts_list = []
        #     for codon in self.codons:
        #         if codon.start <= last and codon.end >= first:
        #             codon_muts = []
                    
        #             for i in xrange(3):     # three bases
        #                 if codon.start + i >= first and codon.start + i <= last: # make sure that mutated position is within the kmer
        #                     codon_muts.extend(list(codon.snv_gen(i)))
                            
        #             codon_muts_list.append(codon_muts)

        #     for muts in itertools.product(*codon_muts_list):
        #         kmer = self.seq[first-self.start:last+1-self.start]
        #         for mut in muts:
        #             kmer = mut.mutate(kmer, first)

        #         yield (kmer, muts)
            
        return
                
        
def group(iterable, key):
    """ group items in an iterable such that the
    dictionary map from the key to the list of items
    """

    d = {}
    for item in iterable:
        registerlist(d, key(item), item)

    return d
    
def group_codons(codons):

    # group condons by chromosomes
    chrm2codons = group(codons, key=lambda x: x.chrm)
    contig_stack = []
    for chrm, subcdns in chrm2codons.iteritems():
        # sort codons by start
        subcdns = sorted(list(subcdns), key=lambda x: x.start)
        
        # push first codon into the contig stack
        contig_stack.append(CodonContig(subcdns[0]))

        # push the rest into the contig stack
        for i in xrange(1, len(subcdns)):
            top = contig_stack[-1]
            if top.end + kmer_len < subcdns[i].start: # non-overlap case
                contig_stack.append(CodonContig(subcdns[i]))
            else:               # overlap case
                top.append(subcdns[i])

    # for chrm in chrm2codons.iteritems():
    # for codon in codons:
    #     print codon, codon.seq, codon.tr_aa(), list(codon.nonsynym_mutations())

    return contig_stack


def retrieve_refseq_by_contigs(referencefile, contigs, extend=0):
    """ every object in contigs must define: chrm, start and end
    this is the memory-intensive version
    """

    for contig in contigs:
        contig._extended_start_ = contig.start - extend
        contig._extended_end_ = contig.end + extend
    
    # group contigs by chromosomes
    chrm2contigs = group(contigs, key=lambda x: x.chrm)
    for chrm, subctgs in chrm2contigs.iteritems():
        # sort contigs by start (reversed)
        chrm2contigs[chrm] = sorted(list(subctgs), key=lambda x: x.start, reverse=True)

    # initialize contig seq
    for contig in contigs:
        contig._seq_ = ''

    with open(referencefile) as f:
        contigs_to_process = []
        for line in f:
            line = line.strip()
            line_length = len(line)
            
            # process chromosome header
            if line.startswith(">"):
                chrm = line[1:]
                curr_contigs = chrm2contigs[chrm] if chrm in chrm2contigs else []
                curr_pos = 1
                assert contigs_to_process == []
                continue

            next_pos = curr_pos + line_length
            
            # process the contigs that start in this line
            while curr_contigs and curr_contigs[-1]._extended_start_ < next_pos:
                contigs_to_process.append(curr_contigs.pop())

            # update the contig seqs
            new_contigs_to_process = []
            for contig in contigs_to_process:
                step = contig._extended_end_ - curr_pos + 1
                if step <= line_length:
                    seg_end = step # contigs that end on this line
                else:
                    seg_end = line_length # contigs that do not end on this line
                    new_contigs_to_process.append(contig)
                    
                seg_start = max(0, contig._extended_start_ - curr_pos)
                contig._seq_ += line[seg_start:seg_end]

            # for contig in contigs_to_process:
            #     if contig not in new_contigs_to_process:
            #         if len(contig._seq_) != 203:
            #             print contig, step
            #             print contig._seq_
            #             print "currpos: ", curr_pos
            #             print "currline: ", line
            #             print "nextpos: ", next_pos
            #             import sys
            #             sys.exit(1)

            contigs_to_process = new_contigs_to_process

            # go to next line
            curr_pos = next_pos

    return

def main():

    # retrieve codon positions
    codons=set()
    with open(os.path.join(pathsetup.sharedir, "mut_nt_pos.txt")) as f:
        for line in f:
            prot, _, chrm, orientation, loc = line.split('\t')

            prot_id, prot_loc = prot.split(":")

            loc = int(loc)

            if orientation == "+":
                codons.add(Codon(chrm, loc, loc+2, '+'))
            else:
                codons.add(Codon(chrm, loc, loc+2, '-'))

    # remove duplicate codons
    codons = list(codons)

    # group codons into contigs
    codon_contigs = group_codons(codons)
    for codon_contig in codon_contigs:
        codon_contig.start -= kmer_len - 1
        codon_contig.end   += kmer_len - 1
    
    print len(codon_contigs)

    # retrieve sequence from hg19.fa
    retrieve_refseq_by_contigs(pathsetup.referenceloc, codons)
    for codon in codons:
        codon.seq = codon._seq_
        print codon, codon.seq, codon.natural_seq(), codon.tr_aa()

    retrieve_refseq_by_contigs(pathsetup.referenceloc, codon_contigs)
    for contig in codon_contigs:
        contig.seq = contig._seq_

    dump(codon_contigs, open(os.path.join(pathsetup.datadir, '2013_10_10_178_codon_contigs.pkl'),'w'))

    return

def main2():
    codon_contigs = load(open(os.path.join(pathsetup.datadir, '2013_10_10_178_codon_contigs.pkl')))

    fout = open('/projects/zchong/data/2013_10_13_kmer_lib','w')
    for i, codon_contig in enumerate(codon_contigs):
        print codon_contig
        print codon_contig.codons
        i = 0
        for kmer, muts in codon_contig.kmer_gen(kmer_len):
            if any([mut.ref != mut.var for mut in muts]):
                fout.write('%s\t%s\n' % (';'.join(['%s:%d:%s->%s(%s->%s)' % (codon_contig.chrm, mut.loc, mut.ref, mut.var, mut.ref_aa, mut.var_aa) for mut in muts]), kmer))

import unittest
class KmerGenTest(unittest.TestCase):

    def testGroupCodon_double_codon(self):
        codon1 = Codon('chr1', 200, 202, '+')
        codon1.seq = 'AAT'
        
        codon2 = Codon('chr1', 207, 209, '+')
        codon2.seq = 'AAT'
        
        codon_contig = group_codons([codon1, codon2])[0]
        self.failUnless(codon_contig.start == 200)
        self.failUnless(codon_contig.end == 209)

    def testGroupCodon_single_codon(self):
        codon1 = Codon('chr1', 200, 202, '+')
        codon1.seq = 'AAT'
        
        codon_contig = group_codons([codon1])[0]
        self.failUnless(codon_contig.start == 200)
        self.failUnless(codon_contig.end == 202)


    def testSNVgen(self):
        codon = Codon('chr1', 200, 202, '+')
        codon.seq = 'AAT'
        codon.snv_gen(1)
        

if __name__ == '__main__':
                         
    # main()
    main2()
    # unittest.main()

# unsorted = loc2codon.keys()
# locs = sorted(unsorted)
# print '\n'.join(map(str, unsorted))

# print '==========='
# print '\n'.join(map(str, locs))

# print loc2codon


        # prev_prot_id = prot_id
        # if var_type == "single base substitution":
        #     locs = map(int, locstr.strip().split(','))
        #     for loc in locs:
        #         if (prot_id, loc) not in prot2muts:
        #             print (prot_id, loc)
        #             raise Exception("protein clear house db unfound")





# for chrm, hotspot in chrm2hotspots.iteritems():
#     print chrm, sorted(hotspot)
