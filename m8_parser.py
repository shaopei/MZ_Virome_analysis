#calculate the number of hit each ref subject has in one sample
#make sure a subject is only counted once in a query.


from sys import argv
from collections import Counter
import glob
import numpy as np
import time


def subject_occurrence_per_file (f):
    # f is file name
    #c is Counter
    with open(f,"U") as data:
        print f
        q_s_pair = Counter()
        query_c = Counter()
        subject_c = Counter()
        for l in data:
            [q,s] = l.strip().split('\t')[0:2]
            if q_s_pair[(q, s)] == 0:
                query_c.update([q])
                subject_c.update([s])
            q_s_pair.update([(q,s)])
    query_hist = Counter(query_c.values())
    return query_hist, subject_c

def subject_occurrence_per_file_quick_memory_demanding (f): #not finished yet
# f is file name
    #data = np.loadtxt(f, dtype=str, delimiter='\t', usecols=xrange(3))
    q_list=[]
    s_list=[]
    for l in open(f,"U").readlines():
        ll=l.strip().split()
        q_list.append(ll[0])
        s_list.append(ll[1])

    query_c = Counter(q_list)
    subject_c = Counter(s_list)
    q_s_pair = Counter(zip(q_list,s_list))

    for qs, count in q_s_pair.most_common():
        if count > 1:
            #print query_c[qs[0]], count
            query_c[qs[0]] = query_c[qs[0]] - count + 1
            subject_c[qs[1]] = subject_c[qs[1]] - count + 1
            #print query_c[qs[0]]
        else:
            #print query_c[qs[0]], count
            break

    query_hist = Counter(query_c.values())
    return query_hist, subject_c





def main():
    start_time = time.time()
    f_query_dic = {}
    f_subject_dic = {}
    for f in glob.glob('*matches.m8'):
        print f
        q, s = subject_occurrence_per_file_quick_memory_demanding (f)
        #q, s = subject_occurrence_per_file(f)
        f_query_dic[f] = q
        f_subject_dic[f] = s
        print 'Elapsed time: %s' %(time.time()-start_time)

    #f_subject_dic = {f:subject_occurrence_per_file (f) for f in glob.glob('*matches.m8')}


    f_list = f_subject_dic.keys()
    f_list.sort()
    query = Counter()
    subject = Counter()
    for f in f_list:
        query = query + f_query_dic[f]
        subject = subject + f_subject_dic[f]

    with open('diamond_blastx_m8_parsed_query.txt',"w") as out:
    #with open('test.txt',"w") as out:
        out.write("\t".join(["query", "total_count"] + [f for f in f_list]))
        out.write("\n")
        l = query.keys()
        l.sort()
        for q in l:
            out.write("\t".join([str(q), str(query[q])] + [str(f_query_dic[f][q]) for f in f_list]))
            out.write("\n")


    with open('diamond_blastx_m8_parsed_subject.txt',"w") as out:
        out.write("\t".join(["subject", "total_count"] + [f for f in f_list]))
        out.write("\n")
        for s, c in subject.most_common():
            out.write("\t".join([s, str(c)] + [str(f_subject_dic[f][s]) for f in f_list]))
            out.write("\n")





if __name__ == '__main__':
    main()