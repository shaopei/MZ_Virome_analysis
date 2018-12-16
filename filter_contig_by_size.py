from sys import argv

print 'python filter_contig_by_size.py cotig.fa size'

f1 = argv[1]
size = argv[2]

seq_header=''
seq_temp=[]
sl=0
#size = 1000

with open(f1+'_LargerThan'+str(size), 'w') as out:
    with open(f1,"U") as f: # automatically close the file when done
        for line in f:#.readlines():
            line = line.strip() # trim endline
            if line.startswith(">"):
                if sl > int(size):
                    out.write(seq_header+'\n')
                    out.write('\n'.join(seq_temp))
                    out.write('\n')
                sl = 0
                seq_header= line
                seq_temp=[]
            else:
                sl += len(line)
                seq_temp.append(line)