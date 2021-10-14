# %%
# try split
import gzip
seq_tar =''
with gzip.open('/data/student2020/resources/hg38.fa.gz') as f:
    for line in f:
        seq_tar = seq_tar + line.decode().strip().replace("\n",'')

#seq_tar = seq_tar.replace("\n", '')
print('fileloaded')
seq_splited = seq_tar.split('chr')
seq_splited = seq_splited[1:]


# read chromosome name
def readname(splitedstr):
    name = []
    if splitedstr.find('alt') == -1:  # 如果没有alt
        count = 0
        Basecount = 0
        ATCG = "ATCGNatcg"
        while (Basecount < 2):
            this = splitedstr[count]
            if ATCG.find(this) != -1:
                Basecount = Basecount + 1;
            else:
                Basecount = 0
            name.append(this)
            count = count + 1
        outcome = 'chr' + ''.join(name[:-2])
        sequence = splitedstr[count - 2:]
    else:
        loca = splitedstr.find('alt') + 3
        outcome = 'chr' + splitedstr[0:loca]
        sequence = splitedstr[loca + 1:]

    return [outcome, sequence]


names = []
seqences = []
for items in seq_splited:
    subsplite = readname(items)
    names.append(subsplite[0])
    seqences.append(subsplite[1])

with open('Name_index.txt', 'w') as q:
    for i in names:
        q.write(i)
        q.write('\n')
with open('overallseq.txt', 'w') as s:
    for i in seqences:
        s.write(i)
        s.write('\n')

print('filesplited')

