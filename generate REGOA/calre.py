import importdatas as datas
from sklearn.metrics.pairwise import cosine_similarity
Function = 'CC'
from multiprocessing import Process, Queue, Manager
import random
import pickle
path = './datamgi/'

def importvec(filename):
    vecs=dict()
    with open(filename,'r') as f:
        next(f)
        for line in f:
            line=line.split()
            vecs[line[0]]=[float(x) for x in line[1:]]
    print('import vecs done')
    return vecs

def getthres():

    fmax = dict()

    terms = set()
    with open(path+Function+'threshold.txt','r') as f:
        next(f)
        for line in f:
            line = line.split()

            fmax[line[0]] = float(line[1])

            terms.add(line[0])
    return fmax, terms

def calres(vecs, res, fmax, goaterms):
    M = len(res)
    count = 0
    prob = dict()

    gfmax = open(path + Function + 're-term-fmax.txt', 'w')

    for re in res:
        prob[re] = dict()
        for term in go.keys():
            prob[re][term] = cosine_similarity([vecs[term], vecs[re]])[0][1]
        terms = set(go.keys())
        while terms:
            term = terms.pop()
            nowmax = prob[re][term]
            nowterm = term
            for tterm in go[term]['descent'] & set(prob[re].keys()):
                if prob[re][tterm] > nowmax:
                    nowmax = prob[re][tterm]
                    nowterm = tterm
            prob[re][term] = nowmax
            for tterm in go[term]['descent'] & go[nowterm]['ancesent'] & terms:
                prob[re][tterm] = nowmax
            terms = terms - go[term]['descent'] & go[nowterm]['ancesent']
        gfmax.write(re+'\t')
        for term in goaterms:
            if prob[re][term] > fmax[term]:
                gfmax.write(term + '\t')

        gfmax.write('\n')

        count += 1
        print(count, '/', M, ' test')

def calress(redict, vecs, res, fmax, goaterms, logo):
    M = len(res)
    count = 0
    prob = dict()
    #g = open('./data' + logo+Function + 're-term-low.txt', 'w')
    for re in res:
        prob[re] = dict()
        #g.write(re+'\t')
        for term in go.keys()&vecs.keys():
            prob[re][term] = cosine_similarity([vecs[term], vecs[re]])[0][1]
        terms = set(go.keys()&vecs.keys())
        while terms:
            term = terms.pop()
            nowmax = prob[re][term]
            nowterm = term
            for tterm in go[term]['descent'] & set(prob[re].keys()):
                if prob[re][tterm] > nowmax:
                    nowmax = prob[re][tterm]
                    nowterm = tterm
            prob[re][term] = nowmax
            for tterm in go[term]['descent'] & go[nowterm]['ancesent'] & terms:
                prob[re][tterm] = nowmax
            terms = terms - go[term]['descent'] & go[nowterm]['ancesent']
        annos = set()
        for term in goaterms:
            if prob[re][term] > fmax[term]:
                annos.add(term)
                #print(term)
                #g.write(term+'\t')
        count += 1
        #g.write('\n')
        redict[re] = annos
        #for term in redict[re]:
            #print(term)
        #print(re + '\t' + '\t'.join(list(redict[re]))+'\n')
        if count % 20 == 0:
            print(logo, count, '/', M, ' test')


def main():
    '''
    global grn
    grn = datas.importgrn()  # return {'genes':{gene:re,gene:re...},'res':{re:{'reg':genes,'bind':tfs}},'tfs':{tf:re,tf:re}}
    print('import grn done')
    '''
    global go
    go = datas.importgo(path+'go.obo')  # {go_id:{'is_a':terms,'children':terms,}}
    go = {x: go[x] for x in go if go[x]['namespace'] == Function}
    vecs = importvec(path+Function+'vecs.txt')
    #res = grn['res'].keys()
    g = open(path + 'REsP3.pkl', 'rb+')
    res = pickle.load(g)
    #res = res[6000:7000]
    fmax, terms = getthres()
    #calres(vecs, res, high, fmax, low, terms)

    redict = Manager().dict()
    N = 8
    resset = []
    for i in range(N):
        resset.append(set())
    for i in res:
        resset[random.randint(0, N - 1)].add(i)
    processes = []
    for i in range(N):
        processes.append(Process(target=calress, args=(redict, vecs, resset[i], fmax, terms, 'p' + str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    with open(path+Function + 're-term-fmax.txt','a') as f:
        for re in redict.keys():
            f.write(re + '\t' + '\t'.join(list(redict[re]))+'\n')

if __name__=='__main__':
    main()