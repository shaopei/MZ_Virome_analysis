from sys import argv
input_otu_table = argv[1]

gene_KO_ann_dic ={}


with open('/workdir/sc2457/IGC_20150727/IGC.annotation_OF.summary', 'U') as fi:
    for l in fi:
        ll = l.strip().split('\t')
        gene_KO_ann_dic[ll[1]] = [ll[11],ll[7]]


ko_ann = open (argv[1][0:-4]+'_ko_ann.txt','w')


f_list = [ko_ann]

def write_to_file (text):
    for f in f_list:
        f.write(text)



with open(input_otu_table, 'U') as f_t:
    for i, l in enumerate(f_t):
        if i == 0:
            write_to_file (l)
        elif i == 1:
            write_to_file (l.strip())
            write_to_file ('\ttaxonomy\n')
        else:
            geneID = l.strip().split('\t')[0]
            write_to_file(l.strip())
            write_to_file('\t')
            try:
                ko_ann.write('; '.join(gene_KO_ann_dic[geneID]))
            except:
                print l
                ko_ann.write('unknown; unknown')
            write_to_file('\n')

for f in f_list:
    f.close()