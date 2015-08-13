# vim: fdm=indent
'''
author:     Fabio Zanini
date:       13/08/15
content:    Make an MSA and tree of haplotypes from all patients, to exclude
            cross-contaminations.
'''
# Modules
import os
import sys
import subprocess as sp 
import argparse
from Bio import SeqIO
from Bio import Phylo

from hivevo.hivevo.patients import Patient
from hivevo.hivevo.HIVreference import HIVreference

from filenames import get_figure_folder



# Functions
def make_tree(region, fn_ali, fn_tree, tmpfile='/tmp/seqs.fasta'):
    '''Make tree of minor haplotype variants from all patients + outgroup'''

    # Collect haplpotypes from patients
    seqs = []
    patients = ['p1', 'p2', 'p3','p5', 'p6', 'p8', 'p9', 'p11']
    for pcode in patients:
        p = Patient.load(pcode)
        
        for seq in p.get_haplotype_alignment(region):
            seq.id = 'patient_'+pcode+'_'+seq.id
            seq.name = 'patient_'+pcode+'_'+seq.name
            seq.description = 'patient '+pcode+', '+seq.description

            seqs.append(seq)

    # Add reference as an outgroup
    ref = HIVreference(load_alignment=False)
    refseq = ref.annotation[region].extract(ref.seq)
    seqs.append(refseq)

    # Align (Muscle)
    if os.path.isfile(tmpfile):
        os.remove(tmpfile)
    SeqIO.write(seqs, tmpfile, 'fasta')

    try:
        sp.call(['muscle', '-maxiters', '1', '-diags', '-in', tmpfile, '-out', fn_ali])
    finally:
        os.remove(tmpfile)

    # Annotate for FastTree (does not accept double labels)
    seqs = []
    for seq in SeqIO.parse(fn_ali, 'fasta'):
        seq.name = seq.name+'_#'+str(len(seqs))
        seq.id = seq.id+'_#'+str(len(seqs))
        seqs.append(seq)
    SeqIO.write(seqs, tmpfile, 'fasta')

    # FastTree
    try:
        sp.call(['FastTree', '-nt', '-out', fn_tree, tmpfile])
    finally:
        os.remove(tmpfile)

    # reroot with outgroup
    tree = Phylo.read(fn_tree, 'newick')
    for leaf in tree.get_terminals():
        if refseq.id in leaf.name:
            break
    tree.root_with_outgroup(leaf)
    Phylo.write(tree, fn_tree, 'newick')





# Script
if __name__=="__main__":

    parser = argparse.ArgumentParser(description="make figure")
    parser.add_argument('--redo', action='store_true',
                        help='recalculate data')
    params=parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'

    regions = ['p17', 'V3']
    for region in regions:
        fn_ali = fn_data + 'haplotype_alignment_crosspatient_'+region+'.fasta'
        fn_tree = fn_data + 'haplotype_tree_crosspatient_'+region+'.newick'

        if (not os.path.isfile(fn_tree)) or params.redo:
            make_tree(region, fn_ali, fn_tree)

        # Check for cross-contamination
        tree = Phylo.read(fn_tree, 'newick')

        # 2. annotate
        def fun_patient(node):
            return node.name.split('_')[1]

        def annotate_tree(node, fun, attrname):
            if not node.clades:
                setattr(node, attrname, [fun(node)])

            else:
                setattr(node, attrname, [])
                for clade in node.clades:
                    annotate_tree(clade, fun, attrname)
                    getattr(node, attrname).extend(getattr(clade, attrname))
                setattr(node, attrname, sorted(set(getattr(node, attrname))))
                    
        annotate_tree(tree.root, fun_patient, 'patient')

        # 3. group by patient
        from collections import defaultdict
        groups = defaultdict(list)
        for node in tree.get_terminals():
            if len(node.patient) == 1:
                groups[node.patient[0]].append(node)

        # 4. check for monophyletic
        print region, 'non monophyletic:',
        non_mp = []
        for pcode, group in groups.iteritems():
            if not tree.is_monophyletic(group):
                non_mp.append(pcode)
        print non_mp

