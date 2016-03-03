import csv

def listInterestMogs():
    interest_mog = {}
    f = open('rprots.txt')
    fr = csv.DictReader(f, delimiter='\t')
    for line in fr:
        interest_mog[line['MOG']] = line['protein'].replace('/','_')
    f.close()
    return interest_mog

def makeMogNamesDict():
    mog_names = {}
    f = open('../info/allMOGShortNames.csv')
    fr = csv.DictReader(f, delimiter='\t')
    for line in fr:
        mog_names[line['MOG']] = line['shortName']
    f.close()
    return mog_names

def findOperons():
# выписываем опероны с интересными MOG
    all_int_operons = {}
    mog_operons = {x:[] for x in interest_mog}
    g = open('../info/allBacteriaGenesInfo2.csv')
    gr = csv.DictReader(g, delimiter='\t')
    for line in gr:
        if line['mogId'] in interest_mog:
            if line['tuId'] not in all_int_operons:
                all_int_operons[line['tuId']] = []
            mog_operons[line['mogId']].append(line['tuId'])
    g.close()

# записываем, какие гены лежат в выписанных оперонах
    g = open('../info/allBacteriaGenesInfo2.csv')
    gr = csv.DictReader(g, delimiter='\t')
    for line in gr:
        if line['tuId'] in all_int_operons:
            all_int_operons[line['tuId']].append(line['mogId'])
    g.close()
    return all_int_operons, mog_operons

def countNeighbs(mog):
# считаем появления каждого гена
    neighbs = {}
    for op in mog_operons[mog]:
        for gene in all_int_operons[op]:
            name = gene+'\t'+interest_mog.get(gene, mog_names[gene])
            neighbs[name] = neighbs.get(name, 0) + 1
    return neighbs

interest_mog = listInterestMogs()
mog_names = makeMogNamesDict()
all_int_operons, mog_operons = findOperons()

for mog in interest_mog:
    neighbs = countNeighbs(mog)
    o = open('neighbs/'+mog+'_'+interest_mog[mog], 'w')
    o.write('MOG\tshortName\tcount\n')
    neighb_tuples = [(x, neighbs[x]) for x in neighbs]
    neighb_tuples.sort(key = lambda x: x[1], reverse=True)
    o.write('\n'.join([x[0]+'\t'+str(x[1]) for x in neighb_tuples]))
    o.close()
