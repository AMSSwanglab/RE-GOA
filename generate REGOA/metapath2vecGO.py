
import importdatas as datas
from random import choice
from gensim.models import Word2Vec
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import precision_recall_curve
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,savefig
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc  ###计算roc和auc
import seaborn as sns
import scipy.cluster.hierarchy as sch
import palettable
from random import shuffle
from multiprocessing import Process, Queue, Manager
import random
import pickle
import calthres

Function = 'CC'
path = './datamgi/'
notename = '-noDis'

def preprocess():



    global grn
    if path[6] == 'h':
        grn = datas.importgrnhg()    #    return {'genes':{gene:re,gene:re...},'res':{re:{'reg':genes,'bind':tfs}},'tfs':{tf:re,tf:re}}
    else:
        grn = datas.importgrnmgi()  # return {'genes':{gene:re,gene:re...},'res':{re:{'reg':genes,'bind':tfs}},'tfs':{tf:re,tf:re}}
    print('import grn done')
    
    global go
    go = datas.importgo(path+'go.obo')      #    {go_id:{'is_a':terms,'children':terms,}}

    go = {x: go[x] for x in go if go[x]['namespace'] == Function}
    print('import go done')
    global goa
    global goaraw
    goaraw = datas.importgoa(go, path+'train.txt')    #    return {'genes':{gene:terms,...},'terms':{term:genes}:}
    print('import goa done')
    goa=datas.goatpr(goaraw, go)
    print('tpr goa done')
    global ppi

    ppi = datas.importppi(path+'ppifull.txt')
    print('import ppi done')

    print('genes num in goa:', len(goa['genes']), 'terms num in goa:', len(goa['terms']))
    #print('genes num in grn:', len(grn['genes']))
    print('terms num in go:', len(go))


    global gene_n_goa
    gene_n_goa = grn['genes'].keys()-goa['genes'].keys()
    print('gene in grn not in goa:', len(gene_n_goa))
    genetf = set(grn['genes'].keys()) | set(grn['tfs'].keys())
    gene_n_grn = goa['genes'].keys()-genetf
    print('gene in goa not in grn:', len(gene_n_grn))

    print('genes num:', len(set(grn['genes'].keys()) & set(goa['genes'].keys()) & set(ppi.keys())))
    print('res num:', len(grn['res']))
    print('tfs num:', len(grn['tfs']))
    print('terms num:', len(go.keys()))


def pickppi(ppisp):
    x=random.randint(0, sum(list(ppisp.values()))-1)
    count = 0
    for item in ppisp.keys():
        count+=ppisp[item]
        if count > x:
            return item
    return item


def genwalk(walks, genesets, walknum, walklen, goa, go, ppi, logo):
    GOAgenes = set(goa['genes'].keys())
    GOAterms = set(goa['terms'].keys())

    M = len(genesets)
    count=0
    for gene in genesets:
        for i in range(walknum):
            walk=[]                                  # GTTTTTTT...
            walk.append(gene)
            term = choice(list(goa['genes'][gene]))
            walk.append(term)
            for j in range(walklen):
                term = choice(list(go[term]['is_a'] | go[term]['children']))

                walk.append(term)
            walks.append(walk)

            walk = []                                #GTGTGTGTGT...
            walk.append(gene)
            term = choice(list(goa['genes'][gene]))
            walk.append(term)
            for j in range(walklen):
                gene = choice(list(goa['terms'][term]))
                walk.append(gene)
                term = choice(list(goa['genes'][gene]))
                walk.append(term)
            walks.append(walk)


            walk = []                                 #GTTGGTTGGTTG...
            walk.append(gene)
            term = choice(list(goa['genes'][gene] & GOAterms))
            walk.append(term)

            tterm = choice(list((go[term]['is_a'] | go[term]['children']) & GOAterms))

            if not goa['terms'][tterm] & GOAgenes & ppi.keys():
                walk.append(tterm)
                walks.append(walk)
                continue
            term = tterm
            walk.append(term)

            tgene = choice(list(goa['terms'][term] & GOAgenes & ppi.keys()))
            if not ppi[tgene].keys() & GOAgenes:
                walk.append(tgene)
                walks.append(walk)
                continue
            gene = tgene
            walk.append(gene)

            for j in range(walklen // 2):
                gene = pickppi({x:ppi[gene][x] for x in ppi[gene].keys() if x in GOAgenes})
                walk.append(gene)
                term = choice(list(goa['genes'][gene] & GOAterms))
                walk.append(term)

                tterm = choice(list((go[term]['is_a'] | go[term]['children']) & GOAterms))
                if not goa['terms'][tterm] & GOAgenes & ppi.keys():
                    walk.append(tterm)
                    break
                term = tterm
                walk.append(term)

                tgene = choice(list(goa['terms'][term] & GOAgenes & ppi.keys()))
                if not ppi[tgene].keys() & GOAgenes:
                    walk.append(tgene)
                    break
                gene = tgene
                walk.append(gene)

                #print(len(walk), walk)
            walks.append(walk)



        count += 1
        if count % 10 == 0:
            print(logo, count, '/goa', M)


def grnwalk(walks, genesets, walknum, walklen, grn, logo):
    M = len(genesets)  # GRTRG
    count = 0
    for gene in genesets:
        for i in range(walknum * 5):
            walk = []
            walk.append(gene)
            for j in range(walklen // 4):
                re = choice(list(grn['genes'][gene]))
                walk.append(re)
                tf = choice(list(grn['res'][re]['binded']))
                walk.append(tf)
                re = choice(list(grn['tfs'][tf]))
                walk.append(re)
                tf = choice(list(grn['res'][re]['reg']))
                walk.append(tf)
            walks.append(walk)
        count += 1
        if count % 10 == 0:
            print(logo, count, 'grn', M)

def mp2vwalk(walknum,walklen):
    global walks
    walks = []
    GOAgenes=set(goa['genes'].keys())
    GOAterms=set(goa['terms'].keys())
    '''
    L=len(grn['res'].keys())          #add re into goa   RE-Term
    count=0
    for re in grn['res'].keys():
        count+=1
        if count%10000==0:
            print(count,'/resgoa',L)
        for gene in (grn['res'][re]['reg']) & GOAgenes:
            for term in goaraw['genes'][gene]:
                walks.append([re,term])

    walks = walks*10
    print(len(walks))
    '''

    g = open(path + Function + 're-goa.pkl', 'rb+')
    res = pickle.load(g)
    for re in res.keys():
        for term in res[re]:
            walks.append([re,term])
    walks = walks * 40
    print(len(walks))
    '''
    tfs = grn['tfs'].keys() & grn['genes'].keys()
    L = len(tfs)
    count = 0
    for tf in tfs:
        count +=1
        print(count, '/ tfs', L)
        for i in range(walknum):
            walk = [tf]
            ttf = tf
            for j in range(walklen):      # TRTRTRT...
                tre = choice(list(grn['genes'][ttf]))
                walk.append(tre)
                if list(grn['res'][tre]['binded'] & tfs):
                    ttf = choice(list(grn['res'][tre]['binded'] & tfs))
                    walk.append(ttf)
                else:
                    break
            walks.append(walk)
    '''
    L=len(grn['res'].keys())
    count=0
    for re in grn['res'].keys():        #RGGRRGGRRGGR...
        count+=1
        if count%100 == 0:
            print(count,'/res',L)
        walk = [re]
        tre = re
        for j in range(walklen//4):
            if not list(grn['res'][tre]['reg'] & ppi.keys()):
                break
            tgene = choice(list(grn['res'][tre]['reg'] & ppi.keys()))
            walk.append(tgene)
            if not list(ppi[tgene].keys() & grn['genes'].keys()):
                break
            tgene = choice(list(ppi[tgene].keys() & grn['genes'].keys()))
            walk.append(tgene)
            tre = choice(list(grn['genes'][tgene]))
            walk.append(tre)
            tre = choice(list(grn['genes'][tgene]))
            walk.append(tre)
        walks.append(walk)
    print(len(walks))

    walksgene = Manager().list()
    N = 8
    geneset = []
    for i in range(N):
        geneset.append(set())
    for i in goa['genes'].keys():
        geneset[random.randint(0, N - 1)].add(i)
    processes = []
    for i in range(N):
        processes.append(Process(target=genwalk, args=(walksgene, geneset[i], walknum, walklen, goa, go, ppi, 'p'+str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    walks = walks+list(walksgene)
    print(len(walks))

    M = len(ppi.keys())              #GGGGGGG (PPI)
    count=0
    for gene in ppi.keys():
        count+=1
        for i in range(walknum):
            walk=[]
            walk.append(gene)
            for j in range(walklen):
                gene=pickppi(ppi[gene])
                walk.append(gene)
            walks.append(walk)
        if count%100==0:
            print(count,'/ppi',M)

    walksgrn = Manager().list()
    N = 8
    geneset = []
    for i in range(N):
        geneset.append(set())
    for i in grn['genes'].keys():
        geneset[random.randint(0, N - 1)].add(i)
    processes = []
    for i in range(N):         # GRTRG
        processes.append(
            Process(target=grnwalk, args=(walksgrn, geneset[i], walknum, walklen, grn, 'p' + str(i))))
    for i in range(N):
        processes[i].start()
    for i in range(N):
        processes[i].join()
    walks = walks + list(walksgrn)
    print(len(walks))
    

    '''
    N=len(grn['genes'].keys())                    #GRTRG
    count=0
    for gene in grn['genes'].keys():
        for i in range(walknum*5):
            walk=[]
            walk.append(gene)
            for j in range(walklen//4):
                re = choice(list(grn['genes'][gene]))
                walk.append(re)
                tf = choice(list(grn['res'][re]['binded']))
                walk.append(tf)
                re = choice(list(grn['tfs'][tf]))
                walk.append(re)
                tf = choice(list(grn['res'][re]['reg']))
                walk.append(tf)
            walks.append(walk)
        count+=1
        if count%100==0:
            print(count,'grn',N)
    '''
    return walks

def main():


    walknum=20
    walklen=50
    #walks=mp2vwalk(walknum,walklen)

    global model
    '''
    f = open('./datahg/diffpaths/' + Function + '_walks' + notename + '.pkl', 'wb')
    pickle.dump(walks, f)
    f.close()
    '''

    g = open('./datahg/diffpaths/' + Function + '_walks' + notename + '.pkl', 'rb+')
    walks = pickle.load(g)
    print('walk done')
    model = Word2Vec(walks, size=128, window=2, min_count=0, sg=1, workers=16, iter=10)
    print('train done')
    model.wv.save_word2vec_format(path+'diffpaths/'+Function+'vecs'+notename+'.txt')


def importvec(filename):
    vecs=dict()
    with open(filename,'r') as f:
        next(f)
        for line in f:
            line=line.split()
            vecs[line[0]]=[float(x) for x in line[1:]]
    print('import vecs done')
    return vecs

def getauroc(y,prob):
    fpr, tpr, threshold = roc_curve(y, prob)  ###计算真正率和假正率
    roc_auc = auc(fpr, tpr)  ###计算auc的值
    lw = 1
    plt.figure(figsize=(10, 10))
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.3f)' % roc_auc)  ###假正率为横坐标，真正率为纵坐标做曲线
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC of metapath2vec in GO&GRN_'+Function)
    plt.legend(loc="lower right")
    plt.savefig(path+"AUROCs/GO&GRN"+Function+".jpg")
    return roc_auc

def test(vecs):
    standard=[]
    score=[]
    with open(path+'test.txt','r') as f:
        count=0
        for line in f:
            line=line.split()
            if line[0] in vecs.keys() and line[1] in vecs.keys()& go.keys():
                #ans=max({model.wv.similarity(line[0],x) for x in go[line[1]]['descent']&set(vecs.keys())}|{model.wv.similarity(line[0],line[1])})
                ans=max({cosine_similarity([vecs[line[0]],vecs[x]])[0][1] for x in go[line[1]]['descent']&set(vecs.keys())}|{cosine_similarity([vecs[line[0]],vecs[line[1]]])[0][1]})
                score.append(ans)
                standard.append(1)
                count+=1
                if count%100==0:
                    print(count,'in 26824 pos')

    with open(path+'neg.txt', 'r') as f:
        count=0
        for line in f:
            line = line.split()
            if line[0] in vecs.keys() and line[1] in vecs.keys()& go.keys():
                #ans=max({model.wv.similarity(line[0],x) for x in go[line[1]]['descent']&set(vecs.keys())}|{model.wv.similarity(line[0],line[1])})
                ans=max({cosine_similarity([vecs[line[0]],vecs[x]])[0][1] for x in go[line[1]]['descent']&set(vecs.keys())}|{cosine_similarity([vecs[line[0]],vecs[line[1]]])[0][1]})
                score.append(ans)
                standard.append(0)
                count+=1
                if count%1000==0:
                    print(count,'in 290000 neg')
    f = open('./datahg/diffpaths/'+Function+'_standard_score'+notename+'.pkl', 'wb')
    pickle.dump((standard,score), f)
    f.close()
    fpr, tpr, threshold = roc_curve(standard,score)  ###计算真正率和假正率
    roc_auc = auc(fpr, tpr)  ###计算auc的值
    print(Function,roc_auc)
    #print(Function, 'AUROC:', getauroc(standard,score))

if __name__=='__main__':
    preprocess()
    main()
    vecs=importvec(path+'diffpaths/'+Function+'vecs'+notename+'.txt')
    test(vecs)
    #calthres.cal(Function, goa, go, vecs, path)