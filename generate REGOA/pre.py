import pickle

Function = 'MF'
dataset = 'fmax'
path = './datahg/'

def importgo(filename=path+'go.obo'):
    # Reading Gene Ontology from OBO Formatted file
    go = dict()
    obj = None
    ns = {'biological_process':'BP','molecular_function':'MF','cellular_component':'CC'}
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
                obj['ancesent']=set()
                obj['descent']=set()
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
                    obj['namespace']=ns[l[1]]
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
            go[go_id]['ancesent']=set()
        temp = list(go[go_id]['is_a'])
        while temp:
            now=temp.pop()
            go[go_id]['ancesent'].add(now)
            temp = list(set(temp)|go[now]['is_a'])
            if 'descent' not in go[now]:
                go[now]['descent']=set()
            go[now]['descent'].add(go_id)
    return go

def main():
    go_back = dict()
    resfmax = dict()
    global go
    go = importgo()
    with open(path+Function+'re-term-'+'fmax'+'.txt','r') as f:
        for line in f:
            line = line.split()
            re = line[0].split('_')
            terms = set(line[1:])
            for term in line[1:]:
                terms = terms | go[term]['ancesent']
            chr = re[0]
            start = int(re[1])
            end = int(re[2])
            if not chr in resfmax.keys():
                resfmax[chr] = dict()
            resfmax[chr][line[0]] = dict()
            resfmax[chr][line[0]]['start'] = start
            resfmax[chr][line[0]]['end'] = end
            resfmax[chr][line[0]]['terms'] = terms
    N = 0
    goanum = 0
    g = open(path+Function+'_REGOA.txt','w')
    res = dict()
    with open(path + Function + 're-term-' + 'raw' + '.txt', 'r') as f:
        for line in f:
            line = line.split()
            re = line[0].split('_')
            terms = set(line[1:])
            for term in line[1:]:
                terms = terms | go[term]['ancesent']
            chr = re[0]
            start = int(re[1])
            end = int(re[2])
            if not chr in res:
                res[chr] = dict()
            if len(terms) == 0:
                continue
            else:
                g.write(line[0]+'\t'+'\t'.join(list(terms & resfmax[chr][line[0]]['terms']))+'\n')
                res[chr][line[0]] = dict()
                res[chr][line[0]]['terms'] = terms & resfmax[chr][line[0]]['terms']
                res[chr][line[0]]['start'] = start
                res[chr][line[0]]['end'] = end
                N += 1
                goanum += len(res[chr][line[0]]['terms'])
                for term in res[chr][line[0]]['terms']:
                    if not term in go_back.keys():
                        go_back[term] = 0
                    go_back[term] += 1
    with open(path+'statistic.txt','a') as f:
        f.write(Function+'\n')
        f.write('RE num:'+ str(N)+'\n')
        f.write('annotation num:'+str(goanum)+'\n')
        f.write('term num:'+str(len(go_back))+'\n')
        f.write('Avg term annotated to per RE:'+str(goanum/N)+'\n')
        f.write('Avg RE annotated by per term:'+str(goanum/len(go_back))+'\n')
    for term in go_back.keys():
        go_back[term] = (go_back[term],N)
        #print(term, go_back[term])
    f = open(path +Function+'re-term-'+'fmaxrawinter'+'.pkl', 'wb')
    pickle.dump(res, f)
    f.close()
    f = open(path + Function + '_backpro_' + 'fmaxrawinter' + '.pkl', 'wb')
    pickle.dump(go_back, f)
    f.close()

def generatethree():
    go_backraw = dict()
    go_backembed = dict()
    go_backgoa = dict()

    resraw = dict()
    resembed = dict()
    resgoa = dict()

    rawgoanum = 0
    embedgoanum = 0
    goanum = 0

    g = open('./datahg/statistic.txt','a')
    g.write(Function + '\n')

    global go
    go = importgo()
    with open(path + Function + 're-term-' + 'fmax' + '.txt', 'r') as f:
        for line in f:
            line = line.split()
            re = line[0].split('_')
            terms = set(line[1:])
            if len(terms) == 0:
                continue
            for term in line[1:]:
                terms = terms | go[term]['ancesent']
            chr = re[0]
            start = int(re[1])
            end = int(re[2])
            if not chr in resembed.keys():
                resembed[chr] = dict()
            resembed[chr][line[0]] = dict()
            resembed[chr][line[0]]['start'] = start
            resembed[chr][line[0]]['end'] = end
            resembed[chr][line[0]]['terms'] = terms
            embedgoanum += len(terms)
            for term in terms:
                if not term in go_backembed.keys():
                    go_backembed[term] = 0
                go_backembed[term] += 1
    REnum = sum([len(chr) for chr in resembed.values()])
    for term in go_backembed.keys():
        go_backembed[term] = (go_backembed[term], REnum)
    g.write('GOAembed\n')

    g.write('REnum:'+str(REnum)+'\n')
    g.write('term num:'+str(len(go_backembed))+'\n')
    g.write('goa num:'+str(embedgoanum)+'\n')
    g.write('Avg term annotated to per RE:' + str(embedgoanum / REnum)+'\n')
    g.write('Avg RE annotated by per term:' + str(embedgoanum / len(go_backembed))+'\n')
    f = open(path + Function + 're-term-' + 'fmax' + '.pkl', 'wb')
    pickle.dump(resembed, f)
    f.close()
    f = open(path + Function + '_backpro_' + 'fmax' + '.pkl', 'wb')
    pickle.dump(go_backembed, f)
    f.close()

    with open(path + Function + 're-term-' + 'raw' + '.txt', 'r') as f:
        for line in f:
            line = line.split()
            re = line[0].split('_')
            terms = set(line[1:])
            if len(terms) == 0:
                continue
            for term in line[1:]:
                terms = terms | go[term]['ancesent']
            chr = re[0]
            start = int(re[1])
            end = int(re[2])
            if not chr in resraw.keys():
                resraw[chr] = dict()
            resraw[chr][line[0]] = dict()
            resraw[chr][line[0]]['start'] = start
            resraw[chr][line[0]]['end'] = end
            resraw[chr][line[0]]['terms'] = terms
            rawgoanum += len(terms)
            for term in terms:
                if not term in go_backraw.keys():
                    go_backraw[term] = 0
                go_backraw[term] += 1
    REnum = sum([len(chr) for chr in resraw.values()])
    for term in go_backraw.keys():
        go_backraw[term] = (go_backraw[term], REnum)
    g.write('GOAraw\n')
    g.write('REnum:'+str(REnum)+'\n')
    g.write('term num:'+str(len(go_backraw))+'\n')
    g.write('goa num:'+str(rawgoanum)+'\n')
    g.write('Avg term annotated to per RE:' + str(rawgoanum / REnum)+'\n')
    g.write('Avg RE annotated by per term:' + str(rawgoanum / len(go_backraw))+'\n')
    f = open(path + Function + 're-term-' + 'raw' + '.pkl', 'wb')
    pickle.dump(resraw, f)
    f.close()
    f = open(path + Function + '_backpro_' + 'raw' + '.pkl', 'wb')
    pickle.dump(go_backraw, f)
    f.close()

    h = open(path + Function + 'REGOA.txt', 'w')
    for chr in resraw.keys():
        resgoa[chr] = dict()
        for re in set(resraw[chr].keys()) & set(resembed[chr].keys()):
            resgoa[chr][re] = dict()
            resgoa[chr][re]['start'] = int(re.split('_')[1])
            resgoa[chr][re]['end'] = int(re.split('_')[2])
            resgoa[chr][re]['terms'] = resraw[chr][re]['terms'] & resembed[chr][re]['terms']
            h.write(re + '\t' + '\t'.join(list(resraw[chr][re]['terms'] & resembed[chr][re]['terms'])) + '\n')
            goanum += len(resraw[chr][re]['terms'] & resembed[chr][re]['terms'])
            for term in resraw[chr][re]['terms'] & resembed[chr][re]['terms']:
                if not term in go_backgoa.keys():
                    go_backgoa[term] = 0
                go_backgoa[term] += 1
    REnum = sum([len(chr) for chr in resgoa.values()])
    for term in go_backgoa.keys():
        go_backgoa[term] = (go_backgoa[term], REnum)
    g.write('GOA\n')
    g.write('REnum:' + str(REnum)+'\n')
    g.write('term num:' + str(len(go_backgoa))+'\n')
    g.write('goa num:' + str(goanum)+'\n')
    g.write('Avg term annotated to per RE:' + str(goanum / REnum)+'\n')
    g.write('Avg RE annotated by per term:' + str(goanum / len(go_backgoa))+'\n')
    f = open(path + Function + 're-term-' + 'fmaxrawinter' + '.pkl', 'wb')
    pickle.dump(resgoa, f)
    f.close()
    f = open(path + Function + '_backpro_' + 'fmaxrawinter' + '.pkl', 'wb')
    pickle.dump(go_backgoa, f)
    f.close()
    g.write('\n\n\n')

if __name__=='__main__':
    #main()
    generatethree()