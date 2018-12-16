import glob
file_names = glob.glob('*_out.fa')
file_names.sort()
with open('Combined_InteMAP_pp_contigs.fa', 'w') as out:
    for fp in file_names:
        prefix = fp.split('_')[0]
        input_fasta = fp
        with open(input_fasta, 'U') as f:
            for l in f:
                if l.startswith('>'):
                    out.write('>'+prefix+'_contig'+l[1:])
                else:
                    out.write(l)