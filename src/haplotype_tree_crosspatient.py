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
import numpy as np
import argparse
from Bio import SeqIO
from Bio import Phylo

from hivevo.hivevo.patients import Patient
from hivevo.hivevo.HIVreference import HIVreference
from util import patient_colors
from filenames import get_figure_folder


# Functions
def make_tree(region, fn_ali, fn_tree, tmpfile='/tmp/seqs.fasta', fasttreebin='FastTree'):
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
        sp.call([fasttreebin, '-nt', '-out', fn_tree, tmpfile])
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
    parser.add_argument('--redo', action='store_true', help='recalculate data')
    parser.add_argument('--fasttreebin', default='FastTree', help='binary of tree builder')
    parser.add_argument('--plot', action='store_true', default=True, help='plot tree')
    params=parser.parse_args()

    username = os.path.split(os.getenv('HOME'))[-1]
    foldername = get_figure_folder(username, 'first')
    fn_data = foldername+'data/'

    regions = ['p17', 'V3', 'IN1']
    for region in regions:
        fn_ali = fn_data + 'haplotype_alignment_crosspatient_'+region+'.fasta'
        fn_tree = fn_data + 'haplotype_tree_crosspatient_'+region+'.newick'

        if (not os.path.isfile(fn_tree)) or params.redo:
            make_tree(region, fn_ali, fn_tree, fasttreebin=params.fasttreebin)

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
                if node.patient[0] in patient_colors:
                    node.color = map(lambda x:int(x*255), patient_colors[node.patient[0]])
                groups[node.patient[0]].append(node)

        # 4. check for monophyletic
        print region, 'non monophyletic:',
        non_mp = []
        for pcode, group in groups.iteritems():
            if not tree.is_monophyletic(group):
                non_mp.append(pcode)
        print non_mp

        # 5. plot the tree colored by patient
        if params.plot:
            import matplotlib.pyplot as plt
            plt.ion()
            for node in tree.get_nonterminals(order='postorder'):
                child_patients = list(set([c.patient[0] for c in node.clades]))
                if len(child_patients)==1 and child_patients[0] in patient_colors:
                    node.color = map(lambda x:int(x*255), patient_colors[child_patients[0]])
                    node.patient = child_patients

            tree.root_at_midpoint()
            tree.ladderize()
            def label_func(x):
                if x.is_terminal() and np.random.random()<0.1:
                    return x.patient[0]
                else:
                    return ''
            fig = plt.figure(figsize = (15,15))
            ax = plt.subplot(111)
            plt.title("Tree of minor variants in all patients from region "+region)
            Phylo.draw(tree, show_confidence=False, label_func=label_func, axes=ax)

            for fmt in ['.png', '.svg', '.pdf']:
                plt.savefig('figures/tree_all_patients_'+region+fmt)
