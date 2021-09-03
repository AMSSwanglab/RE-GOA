def importgoa(go,filename):
    genes=dict()
    terms=dict()
    with open (filename,'r') as f:
        for line in f:
            line = line.split()
            gene = line[0]
            term = line[1]
            if term in go:
                if not gene in genes:
                   genes[gene]=set()
                genes[gene].add(term)
                if not term in terms:
                    terms[term]=set()
                terms[term].add(gene)

    return {'genes':genes,'terms':terms}#,'genetpr':genetpr,'termtpr':termtpr

def goatpr(goa,go): #true path rule
    genes=goa['genes']
    terms=goa['terms']
    genetpr=dict()
    termtpr=dict()
    for term in terms.keys():
        termtpr[term] = terms[term]
    for gene in genes.keys():
        genetpr[gene]=genes[gene]
        for term in list(genes[gene]):
            for t in go[term]['ancesent']:
                genetpr[gene].add(t)
                if not t in termtpr:
                    termtpr[t]=set()
                termtpr[t].add(gene)
    return {'genes':genetpr,'terms':termtpr}#  'genetpr':genetpr,'termtpr':termtpr

def importgrnhg(filename='./datahg/grn-human.txt'):
    genes=dict()
    res=dict()
    tfs=dict()
    with open(filename , 'r') as f:
        next(f)
        for line in f:
            line=line.split()
            gene=line[0]
            re=line[1]
            tf=line[2]
            if not gene in genes:
                genes[gene]=set()
            genes[gene].add(re)
            if not re in res:
                res[re]={'reg':set(),'binded':set()}
            res[re]['reg'].add(gene)
            res[re]['binded'].add(tf)
            if not tf in tfs:
                tfs[tf] = set()
            tfs[tf].add(re)
    return {'genes':genes,'res':res,'tfs':tfs}

def importgrnmgi(filename='./datamgi/grn-mouse.txt'):
    genes=dict()
    res=dict()
    tfs=dict()
    with open(filename,'r') as f:
        for line in f:
            line=line.split()
            gene=line[0]
            re=line[2]
            tf=line[-2].split(';')
            if not gene in genes:
                genes[gene]=set()
            genes[gene].add(re)
            if not re in res:
                res[re]={'reg':set(),'binded':set()}
            res[re]['reg'].add(gene)
            for tff in tf:
                res[re]['binded'].add(tff)
                if not tff in tfs:
                    tfs[tff]=set()
                tfs[tff].add(re)
    return {'genes':genes,'res':res,'tfs':tfs}

def importgo(filename):
    # Reading Gene Ontology from OBO Formatted file
    go = dict()
    obj = None
    ns={'biological_process':'BP','molecular_function':'MF','cellular_component':'CC'}
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

def importppi(filename):
    s=dict()
    with open(filename,'r') as f:
        for line in f:
            line=line.split()
            if not line[0] in s:
                s[line[0]]=dict()
            if not line[1] in s:
                s[line[1]]=dict()
            weight = int(line[2])
            weight = 1
            s[line[0]][line[1]]=weight
            s[line[1]][line[0]]=weight
    return s