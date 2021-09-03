from metapath2vecGO import importvec
from importdatas import importgo
from sklearn.metrics.pairwise import cosine_similarity




def test(vecs):
    haven = set(vecs.keys())
    target = dict()
    functions = go.keys()
    func_index = dict()
    for i, go_id in functions:
        func_index[go_id] = i
    functions = set(functions)
    with open('./data/test.txt','r') as f:
        for line in f:
            line = line.split()
            if line[0] in haven and line[1] in haven:
                if line[0] not in target:
                    target[line[0]] = set()
                target[line[0]].add(line[1])
    labels = list()
    preds = list()
    print(len(target.keys()))
    for gene in target.keys():
        label = [0] * len(functions)
        pred = [0] * len(functions)
        for term in functions:
            pred[func_index[term]] = max({cosine_similarity([vecs[gene],vecs[x]])[0][1] for x in go[term]['descent']&set(vecs.keys())}|{cosine_similarity([vecs[gene],vecs[term]])[0][1]})
        for term in target[gene]:
            label[term] = 1
            for ans in go[term]['anscent']:
                label[ans] = 1



def main():
    vecs = importvec('./data/vecs.txt')
    global go
    go = importgo()  # {go_id:{'is_a':terms,'children':terms,}}
    go = {x: go[x] for x in go if go[x]['namespace'] == 'BP'}
    print('import go done')

    test(vecs)


if __name__=='__main__':
    main()