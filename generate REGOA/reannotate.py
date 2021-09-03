import importdatas as datas
import pickle
Functions = {'MF', 'BP', 'CC'}

def main():
    global go
    go = datas.importgo('./datahg/go.obo')
    global goa
    global goaraw
    goaraw = datas.importgoa(go, './datahg/full.txt')  # return {'genes':{gene:terms,...},'terms':{term:genes}:}
    print('import goa done')
    goa = datas.goatpr(goaraw, go)
    print('tpr goa done')
    global grn
    grn = datas.importgrnhg()  # return {'genes':{gene:re,gene:re...},'res':{re:{'reg':genes,'bind':tfs}},'tfs':{tf:re,tf:re}}
    print('import grn done')
    res = list(grn['res'].keys())
    
    f = open('./datahg/REsP1.pkl', 'wb')
    pickle.dump(res[:50000], f)
    f.close()
    f = open('./datahg/REsP2.pkl', 'wb')
    pickle.dump(res[50000:100000], f)
    f.close()
    f = open('./datahg/REsP3.pkl', 'wb')
    pickle.dump(res[100000:], f)
    f.close()
    

    gopmf = {x for x in go.keys() if go[x]['namespace'] == 'MF'}
    gopbp = {x for x in go.keys() if go[x]['namespace'] == 'BP'}
    gopcc = {x for x in go.keys() if go[x]['namespace'] == 'CC'}

    outmf = open('./datahg/MFre-term-raw.txt', 'w')
    outbp = open('./datahg/BPre-term-raw.txt', 'w')
    outcc = open('./datahg/CCre-term-raw.txt', 'w')

    rerawbp = dict()
    rerawmf = dict()
    rerawcc = dict()
    L = len(grn['res'].keys())  # add re into goa   RE-Term
    count = 0
    notermre = 0
    for re in grn['res'].keys():
        count += 1
        if count % 10000 == 0:
            print(count, '/resgoa', L)
        terms = set()
        if not grn['res'][re]['reg'] & goa['genes'].keys():
            notermre +=1

        for gene in grn['res'][re]['reg'] & goa['genes'].keys():
            terms = terms | goaraw['genes'][gene]
        outmf.write(re)
        outbp.write(re)
        outcc.write(re)
        rerawbp[re] = terms & gopbp
        rerawmf[re] = terms & gopmf
        rerawcc[re] = terms & gopcc

        for term in terms & gopmf:
            outmf.write('\t'+term)
        for term in terms & gopbp:
            outbp.write('\t'+term)
        for term in terms & gopcc:
            outcc.write('\t'+term)
        outmf.write('\n')
        outbp.write('\n')
        outcc.write('\n')
        
    f = open('./datahg/BPre-goa.pkl', 'wb')
    pickle.dump(rerawbp, f)
    f.close()
    f = open('./datahg/MFre-goa.pkl', 'wb')
    pickle.dump(rerawmf, f)
    f.close()
    f = open('./datahg/CCre-goa.pkl', 'wb')
    pickle.dump(rerawcc, f)
    f.close()
    print(notermre)



if __name__=='__main__':
    main()