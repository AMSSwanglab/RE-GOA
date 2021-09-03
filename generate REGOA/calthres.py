from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd


def testterm(prob,stand):
    threshigh, thresfmax, threslow = 1,1,1
    fmax = 0
    findhigh = False
    findlow = False
    p = sum(stand)
    if p==0:
        return 1,1,1
    for t in range(1,1000):
        thres  = t/1000
        tp = 0
        fp = 0
        for x in range(len(prob)):
            if prob[x]>thres:
                if stand[x] == 1:
                    tp+=1
                else:
                    fp+=1
        rec = tp / p
        if tp+fp != 0 :
            pre = tp/(tp+fp)
            if pre >0.3 and not findlow:
                findlow = True
                threslow = thres
            if pre >0.7 and not findhigh:
                findhigh = True
                threshigh = thres
            if pre + rec != 0 :
                f1 = 2*pre*rec/(pre+rec)
                if f1 > fmax:
                    fmax = f1
                    thresfmax = thres
    return threshigh,thresfmax,threslow

def cal(function, goa, go, vecs):
    global functions
    goterm = list(goa['terms'].keys())
    prob = dict()
    standard = dict()

    proteins = set()
    with open('./data/validpro.txt','r') as f:
        count = 0
        for line in f:
            line = line.split()[0]
            if line in vecs.keys():
                prob[line] = dict()
                standard[line] = dict()
                proteins.add(line)
            count+=1


    M = len(prob)
    count = 0
    for protein in prob:
        for term in goterm:
            prob[protein][term] = cosine_similarity([vecs[term] , vecs[protein]])[0][1]
            standard[protein][term] = 0
        terms = set(go.keys())
        while terms:
            term = terms.pop()
            nowmax = prob[protein][term]
            nowterm = term
            for tterm in go[term]['descent'] & set(prob[protein].keys()):
                if prob[protein][tterm] > nowmax:
                    nowmax = prob[protein][tterm]
                    nowterm = tterm
            prob[protein][term] = nowmax
            for tterm in go[term]['descent'] & go[nowterm]['ancesent'] & terms :
                prob[protein][tterm] = nowmax
            terms = terms - go[term]['descent'] & go[nowterm]['ancesent']
        count+=1
        print(count,'/',M,' test')

    with open('./data/valid.txt','r') as f:
        count = 0
        for line in f:
            line = line.split()
            if line[0] in prob and line[1] in goterm:
                standard[line[0]][line[1]] = 1
                for term in set(goterm)&go[line[1]]['ancesent']:
                    standard[line[0]][term] = 1
            count +=1
            print(count,'/goa test')

    with open('./data/'+function+'threshold.txt','w') as f:
        f.write('go_id\t pre=0.7\t Fmax\t pre=0.3\n')
        count = 0
        M=len(goterm)
        for term in goterm:
            threshigh ,thresfmax,threslow = testterm([prob[protein][term] for protein in proteins], [standard[protein][term] for protein in proteins])
            f.write(term+'\t'+str(threshigh)+'\t'+str(thresfmax)+'\t'+str(threslow)+'\n')
            count+=1
            print(count,'/',M,'term test')
