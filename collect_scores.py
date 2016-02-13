# -*- coding: utf-8 -*-
# скрипт берёт: файл-список интересующих MOG; файл про все гены, файл с координатами апстримов, файлы со скорами
# скрипт выдаёт для каждого MOG: файл со скором в кажом геноме
import os.path

# делаем словарь с интересующими нас MOG и их именами
mog_names = {}
f = open('rgenes/rprots.txt')
for line in f:
    p = line.rstrip('\n').split()
    mog_names[p[1]] = p[0]
f.close()
print 'mog_names', len(mog_names)

# делаем словарь с координатами апстримов для всех генов
coords = {}
h = open('info/upstreams-350+50-80')
h.readline()
for line in h:
    p = line.rstrip('\n').split()
    coords[p[0]] = [int(p[2]), int(p[3]), p[1]]
h.close()
print 'coords', len(coords)

# выписываем номера оперонов, в которых лежат нужные MOG, и номера таксонов
taxids = []
operons = []
g = open('info/allBacteriaGenesInfo2.csv')
for line in g:
    p = line.rstrip('\n').split()
    tax, mog, op = p[1], p[5], p[6]
    if mog in mog_names:
        if op not in operons:
            operons.append(op)
        if tax not in taxids:
            taxids.append(tax)
g.close()
print 'taxids', len(taxids)
print 'operons', len(operons)

# создаём словарь геном:список_координат всех генов из оперонов, попавших в operons
# такие списки можно будет отсортировать, чтобы искать скоры
a = {x:[] for x in taxids}
g = open('info/allBacteriaGenesInfo2.csv')
for i in g:
    p = i.rstrip('\n').split()
    gene, tax, mog, op = p[0], p[1], p[5], p[6]
    # if op in op_mogs and tax == '511145':
    if op in operons:
        if gene in coords:
            c = coords[gene]+[gene, mog, op]
            a[tax].append(c)
g.close()
print 'a',len(a)

# сортируем по первой координате
for tax in a:
    a[tax].sort(key=lambda x: x[0])

# извлекаем скоры по координатам
strand = {'+':'false', '-':'true'}
unavailable = []
for tax in a:
    fname = '/mnt/mapr/user/chezoya/genomeScores/'+tax
    if os.path.isfile(fname):
        s = open(fname)
        # s = open(tax)
        n = 0
        scores = []
        for line in s:
            if n < len(a[tax]):
                p = line.rstrip('\n').split()
                if strand[p[5]] == a[tax][n][2] and int(p[1]) >= a[tax][n][0] and int(p[2]) <= a[tax][n][1]:
                    scores.append(float(p[4]))
                if int(p[2]) > a[tax][n][1]:
                    if scores:
                        a[tax][n].append(min(scores))
                    else:
                        a[tax][n].append(0)
                    print n, scores
                    scores = []
                    n += 1
            else:
                break
        s.close()
        if n < len(a[tax]):
            if scores:
                a[tax][n].append(min(scores))
                scores = []
                n += 1
            for i in range(n, len(a[tax])):
                a[tax][i].append(0)
    else:
        print 'scores for', tax, 'are unavailable'
	unavailable.append(tax)

for i in unavailable:
    del a[i]

# теперь нет необходимости держать списки про каждый из генов, сделаем словарь
# запишем все найденные скоры в файлы по taxID (чтобы_было)
for tax in a:
    o = open('rgenes/interScores/'+tax, 'a')
    o.write('locusID\tMOG\ttuID\tscore\n')
    for gene in a[tax]:
        if len(gene)!=7:
            print 'tax', tax
            print 'gene', gene
            print a[tax]
        o.write(gene[3]+'\t'+gene[4]+'\t'+gene[5]+'\t'+str(gene[6])+'\n')
    a[tax] = {x[3]:[x[4],x[6]] for x in a[tax]}

# создаём словарь предыдущих в опреоне генов
pio = {}
k = open('info/PIO_all')
for line in k:
    p = line.rstrip('\n').split()
    if p[0] in pio:
        pio[p[0]].append(p[1])
    else:
        pio[p[0]] = [p[1]]
k.close()

# выбираем из скоров в опероне только те, что лежат выше гена, минимальный записываем в файл
for tax in a:
    for gene in a[tax]:
        mog = a[tax][gene][0]
        if mog in mog_names:
            s = []
            for i in pio[gene]:
                if i in a[tax]:
                    s.append(a[tax][gene][1])
            o = open('rgenes/operonScores/'+mog, 'a')
            o.write(tax+'\t'+str(min(s))+'\n')
            o.close()
