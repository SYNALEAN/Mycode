import pandas as pd
import copy
import json
import numpy as np

with open("Name_index.txt", "r") as f:  # 打开文件
    names = f.readlines()

for line in range(0, len(names)):
    names[line] = names[line].replace("\n", '')  # 去掉列表中每一个元素的换行符 #去掉列表中每一个元素的换行符
print(str(len(names)) + 'names; finished reading')

inputstr = ''
with open("/data/student2020/resources/HBB_genome.fasta", "r") as f:
    next(f)
    for line in f:
        inputstr = inputstr + line.replace("\n", '')

# with open('11seq.txt','r') as f:
#    for line in f:
#        inputstr = inputstr + line.replace("\n", '')

inputstr = inputstr.replace("\n", '')
inputstr = inputstr.replace('A', 't')
inputstr = inputstr.replace('C', 'g')
inputstr = inputstr.replace('G', 'C')
inputstr = inputstr.replace('T', 'A')
inputstr = inputstr.replace('g', 'G')
inputstr = inputstr.replace('t', 'T')

inputstr = inputstr[::-1]


def getMatching(Thefrag, refer):  # 要在这一步加filter吗 #append 太慢了  听所的建议是创建列表再最后转成df
    tarLoc = []
    genoLoc = []

    for i in range(0, len(Thefrag) - 11):
        window_ = Thefrag[i:i + 11]
        # print("testing")
        if window_ in refer:
            #    print("running")
 #           print(window_)
            toadd = copy.copy(refer[window_])
            for item in toadd:
                tarLoc.append(i)
                genoLoc.append(item)
                #    toadd.append(i)
        #    output.append(toadd)
    output = pd.DataFrame({'tarLoc': tarLoc, 'genoLoc': genoLoc})
    return output


def filter(df):
    outputmix = np.zeros((10, 4))
    thedf = df
    b = thedf['genoLoc'] - thedf['tarLoc']
    b = b.astype('int')
    b = b / 1000
    b = b.astype('int')
    thedf['dis'] = b
    cc = b.value_counts()
#    print("pinlv")
#    print(cc)
    index = cc.index
    for i in range(0, min(10,len(index))):#if not enough!modified
        peak = index[i]
        chosed = df[df['dis'] == peak]
        maxrange = chosed.max()
        minrange = chosed.min()
        outputmix[i, 0] = minrange.iloc[0]  # tarloc
        outputmix[i, 1] = maxrange.iloc[0]
        outputmix[i, 2] = minrange.iloc[1]  # genoloc
        outputmix[i, 3] = maxrange.iloc[1]

    outputmix = outputmix.astype('int')
    return outputmix


def matchScore(cha1, cha2, points):
    a = 0
    if cha1 == cha2:
        a = points
    else:
        a = -2
    return a


def SW_match(seq1, seq2):
    scoreMatrix = np.zeros((len(seq1) + 1, len(seq2) + 1))
    # because the initial
    # for i in range(0,len(seq1)+1):
    #     scoreMatrix[i,0]=-2*i
    # for i in range(0,len(seq2)+1):
    #     scoreMatrix[0,i]=-2*i
    # build the matrix,
    G_maxScroe = 0
    LocY = 0
    LocX = 0
    for m in range(1, len(seq1) + 1):
        for n in range(1, len(seq2) + 1):
            mScore = matchScore(seq1[m - 1], seq2[n - 1], 4)
            maxscore = max(scoreMatrix[m - 1, n - 1] + mScore, scoreMatrix[m - 1, n] - 2, scoreMatrix[m, n - 1] - 2,
                           0)  # panishment and 0
            scoreMatrix[m, n] = maxscore
            if maxscore > G_maxScroe:
                G_maxScroe = maxscore
                LocY = m  # 这个是谁的
                LocX = n
            # back trace
    stepy = LocY
    stepx = LocX
    alied_seq1 = []
    alied_seq2 = []
    while scoreMatrix[stepy, stepx] != 0:
        backone = max(scoreMatrix[stepy, stepx - 1], scoreMatrix[stepy - 1, stepx - 1], scoreMatrix[stepy - 1, stepx])
        if scoreMatrix[stepy - 1, stepx - 1] == backone:
            alied_seq1.insert(0, seq1[stepy - 1])
            alied_seq2.insert(0, seq2[stepx - 1])
            stepy = stepy - 1
            stepx = stepx - 1
        elif scoreMatrix[stepy, stepx - 1] == backone:
            alied_seq1.insert(0, '-')
            alied_seq2.insert(0, seq2[stepx - 1])
            stepx = stepx - 1
        else:
            alied_seq1.insert(0, seq1[stepy - 1])
            alied_seq2.insert(0, '-')
            stepy = stepy - 1

    score_per_bite = G_maxScroe / len(alied_seq1)
    return [alied_seq1, alied_seq2, G_maxScroe, stepy, stepx]


def steamlineChr(chro, inputstr, names):
    # 先只看染色体11  #现在要写到全部里去了哦
    print("start " + names[chro])
    with open(str(chro) + '.txt', "r") as f:
        sequence = f.read()
    sequence = sequence.upper()

    with open('./json/' + str(chro) + '_reference.json', 'r') as f:
        reference = json.load(f)

    #    print("All DATA LOADED")

    loclist = getMatching(inputstr, reference)

    #    print(loclist)

    chosedregion = filter(loclist)
    print(chosedregion)

    canlist = []

    # mattching and

    enlarger = 20
    for index in range(0, 10):
        tarstar = max(chosedregion[index, 0] - enlarger, 0)
        tarend = min(chosedregion[index, 1] + enlarger, len(inputstr))
        genostar = max(chosedregion[index, 2] - enlarger, 0)
        genoend = min(chosedregion[index, 3] + enlarger, len(sequence))
        targetaria = inputstr[tarstar:tarend]
        genomearia = sequence[genostar:genoend]
        # 算分吧
        candidate = SW_match(targetaria, genomearia)
        candidate.append([names[chro], tarstar + candidate[3], genostar + candidate[4]])
        del candidate[3:5]
        candidate[0] = ''.join(candidate[0])
        candidate[1] = ''.join(candidate[1])
        # 染色体-目标序列的位置，在基因组的位置
        canlist.append(candidate)

    print("Preparing for printing")
    canlist = pd.DataFrame(canlist, columns=['tar', 'gen', 'score', 'chrLoc'])
    canlist = canlist.sort_values(by="score", ascending=False)
    canlist['tar'] = canlist['tar'].astype(str)
    canlist['gen'] = canlist['gen'].astype(str)

    return canlist


def formatop(canlist):
    output = []
    for i in range(0, 10):
        infolist = canlist.iloc[i, 3]
        strings_ = ''
        tarseq = canlist.iloc[i, 0]
        genseq = canlist.iloc[i, 1]
        middle = ''
        pre_l = 0
        countmatch = 0
        for l in range(0, len(tarseq)):
            if tarseq[l] == genseq[l]:
                middle = middle + '|'
                countmatch = countmatch+1
            else:
                middle = middle + ' '
            if l % 79 == 0 or l == len(tarseq) - 1:
                strings_ = strings_ + '\n'
                strings_ = strings_ + tarseq[pre_l:l] + '\n'
                strings_ = strings_ + middle[pre_l:l] + '\n'
                strings_ = strings_ + genseq[pre_l:l] + '\n'

                pre_l = l
        info_ = '\n' + infolist[0] + '\t' + 'Target Starts at: ' + str(infolist[1]) + \
                   '\t' + 'Mapped to: ' + str(infolist[2]) +'\tMatch: '+ str(countmatch)+'\\'+\
                str(len(tarseq))+'\tScore: ' + str(canlist.iloc[i, 2])

        strings_ = info_ + strings_
        print(strings_)
        # output.append(strings_)


# 对每个染色体都做上面的事情，然后选score的top10


if __name__ == '__main__':

    import multiprocessing

    chrolist = names
    presult=[]
    result = pd.DataFrame(columns=['tar', 'gen', 'score', 'chrLoc'])
    print('pool')
    print(names)
    pool = multiprocessing.Pool(6)
    for i in range(len(chrolist)):
        res= pool.apply_async(steamlineChr,(i,inputstr,names,))
        presult.append(res)
    pool.close()
    pool.join()
    print(presult[0].get())
    for item in presult:
      #  print("running")
        result=result.append(item.get())
    result = result.sort_values(by="score", ascending=False)
    result = result.reset_index(drop=True)
    result = result.iloc[:10,]
    print(result)



    formatop(result)
