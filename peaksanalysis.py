from scipy.stats import binom_test
import pickle
from math import log10


def calculate(goset, go_back, N, resultpath):
    gosetp = dict()
    for term in goset.keys():
        if goset[term] / N > go_back[term][0] / go_back[term][1]:
            gosetp[term] = binom_test(goset[term], N, go_back[term][0] / go_back[term][1])

    terms = list(gosetp.keys())
    terms.sort(key=lambda x: gosetp[x])
    i = 0
    gosetq = dict()
    for term in terms:
        i += 1
        p = gosetp[term]
        gosetq[term] = p * len(terms) / i

    terms.sort(key=lambda x: gosetq[x])
    g = open(resultpath, 'w')
    for term in terms:
        if gosetp[term] > 0.05:
            break
        g.write(term + '\t' + go[term]['name'])
        g.write('\t' + str(goset[term]) + '\t' + str(N) + '\t' + str(go_back[term][0] / go_back[term][1]))
        g.write('\t' + str(gosetp[term]) + '\t')
        if gosetp[term] != 0:
            g.write(str(-log10(gosetp[term])))
        g.write('\t' + str(gosetq[term]) + '\t')
        if gosetq[term] != 0:
            g.write(str(-log10(gosetq[term])))
        g.write('\n')
    g.close()
    return 0


def main(bedfile, dataset, resultpath):
    g = open(dataset + 'MF' + 're-term-fmaxrawinter.pkl', 'rb+')
    resmf = pickle.load(g)
    g = open(dataset + 'MF' + '_backpro_fmaxrawinter.pkl', 'rb+')
    go_backmf = pickle.load(g)

    g = open(dataset + 'BP' + 're-term-fmaxrawinter.pkl', 'rb+')
    resbp = pickle.load(g)
    g = open(dataset + 'BP' + '_backpro_fmaxrawinter.pkl', 'rb+')
    go_backbp = pickle.load(g)

    g = open(dataset + 'CC' + 're-term-fmaxrawinter.pkl', 'rb+')
    rescc = pickle.load(g)
    g = open(dataset + 'CC' + '_backpro_fmaxrawinter.pkl', 'rb+')
    go_backcc = pickle.load(g)

    g = open(dataset + 'GRN.pkl', 'rb+')
    grn = pickle.load(g)

    geneset = dict()

    N = 0
    res = dict()
    for chr in resmf:
        res[chr] = dict()
        for re in resmf[chr]:
            res[chr][re] = resmf[chr][re]
        for re in resbp[chr]:
            res[chr][re] = resbp[chr][re]
        for re in rescc[chr]:
            res[chr][re] = rescc[chr][re]

    gosetmf = dict()
    gosetbp = dict()
    gosetcc = dict()
    with open(bedfile, 'r') as f:
        for line in f:
            line = line.split()
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            dismin = 99999999
            if not chr in res:
                continue
            for re in res[chr].keys():
                if res[chr][re]['end'] < start:
                    dis = start - res[chr][re]['end']
                elif res[chr][re]['start'] > end:
                    dis = res[chr][re]['start'] - end
                else:
                    dismin = 0
                    cloest_re = re
                    break
                if dis < dismin:
                    dismin = dis
                    cloest_re = re
            distances.append(dismin)
            if dismin < 1000:
                N += 1
                if cloest_re in resmf[chr]:
                    for term in resmf[chr][cloest_re]['terms']:
                        if not term in gosetmf:
                            gosetmf[term] = 0
                        gosetmf[term] += 1
                if cloest_re in resbp[chr]:
                    for term in resbp[chr][cloest_re]['terms']:
                        if not term in gosetbp:
                            gosetbp[term] = 0
                        gosetbp[term] += 1
                if cloest_re in rescc[chr]:
                    for term in rescc[chr][cloest_re]['terms']:
                        if not term in gosetcc:
                            gosetcc[term] = 0
                        gosetcc[term] += 1
                for gene in grn['res'][cloest_re]['reg']:
                    if not gene in geneset:
                        geneset[gene] = 0
                    geneset[gene] += 1
    calculate(gosetmf, go_backmf, N, resultpath + '_MF.txt')
    calculate(gosetbp, go_backbp, N, resultpath + '_BP.txt')
    calculate(gosetcc, go_backcc, N, resultpath + '_CC.txt')

    genes = list(geneset.keys())
    genes.sort(key=lambda x: geneset[x], reverse=True)
    with open(resultpath + '_genes.txt', 'w') as f:
        for gene in genes:
            f.write(gene + '\t' + str(geneset[gene]) + '\n')
    return 0


def importgo(filename='./go.obo'):
    # Reading Gene Ontology from OBO Formatted file
    go = dict()
    obj = None
    ns = {'biological_process': 'BP', 'molecular_function': 'MF', 'cellular_component': 'CC'}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == '[Term]':
                if obj is not None:
                    go[obj['id']] = obj
                obj = dict()
                obj['is_a'] = set()
                obj['part_of'] = set()
                obj['regulates'] = set()
                obj['is_obsolete'] = False
                obj['ancesent'] = set()
                obj['descent'] = set()
                continue
            elif line == '[Typedef]':
                if obj is not None:
                    go[obj['id']] = obj
                obj = None
            else:
                if obj is None:
                    continue
                l = line.split(": ")
                if l[0] == 'id':
                    obj['id'] = l[1]
                elif l[0] == 'is_a':
                    obj['is_a'].add(l[1].split(' ! ')[0])
                elif l[0] == 'name':
                    obj['name'] = l[1]
                elif l[0] == 'is_obsolete' and l[1] == 'true':
                    obj['is_obsolete'] = True
                elif l[0] == 'namespace':
                    obj['namespace'] = ns[l[1]]
    if obj is not None:
        go[obj['id']] = obj
    for go_id in list(go.keys()):
        if go[go_id]['is_obsolete']:
            del go[go_id]
    for go_id, val in go.items():
        if 'children' not in val:
            val['children'] = set()
        for p_id in val['is_a']:
            if p_id in go:
                if 'children' not in go[p_id]:
                    go[p_id]['children'] = set()
                go[p_id]['children'].add(go_id)
    for go_id in go.keys():
        if 'ancesent' not in go[go_id]:
            go[go_id]['ancesent'] = set()
        temp = list(go[go_id]['is_a'])
        while temp:
            now = temp.pop()
            go[go_id]['ancesent'].add(now)
            temp = list(set(temp) | go[now]['is_a'])
            if 'descent' not in go[now]:
                go[now]['descent'] = set()
            go[now]['descent'].add(go_id)
    return go


if __name__ == '__main__':
    global go
    go = importgo()
    path = './inputfiledict/'
    datapath = './datamgi/'
    resultpath = './outputfiledict/'
    count = 0
    file = 'exampleinput.bed'

    with open('./configures.txt','r') as f:
        sets = list()
        for line in f:
            sets.append(line.split()[1])
        path = sets[0]
        datapath = sets[1]
        resultpath = sets[2]
        file = sets[3]
    global distances
    distances = list()

    main(bedfile=path + file, dataset=datapath, resultpath=resultpath + file[:-4])
