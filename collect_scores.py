# -*- coding: cp1251 -*-
# скрипт берёт: файл-список интересующих MOG; файл про все гены, файл с координатами апстримов, файлы со скорами
# скрипт выдаёт для каждого MOG: файл со скором в кажом геноме

# делаем словарь с интересующими нас MOG и их именами
mog_names = {}
f = open('ecoliRprotMOG.txt')
for line in f:
    p = line.rstrip('\n').split('\t')
    mog_names[p[1]] = [p[0]]
f.close()

# делаем словарь с координатами апстримов для всех генов
coords = {}
h = open('../info/upstreams-350+50-80')
h.readline()
for line in h:
    p = line.rstrip('\n').split()
    coords[p[0]] = [int(p[2]), int(p[3]), p[1]]
h.close()

# выписываем номера оперонов, в которых лежат нужные MOG, и номера таксонов
taxids = []
operons = []
g = open('../info/allBacteriaGenesInfo2.csv')
for line in g:
    p = line.rstrip('\n').split()
    tax, mog, op = p[1], p[5], p[6]
    if mog in mog_names and tax == '511145':
        if op not in operons:
            operons.append(op)
        if tax not in taxids:
            taxids.append(tax)
g.close()

# создаём словарь геном:список_координат всех генов из оперонов, попавших в operons
# такие списки можно будет отсортировать, чтобы искать скоры
a = {x:[] for x in taxids}
g = open('../info/allBacteriaGenesInfo2.csv')
for i in g:
    p = i.rstrip('\n').split()
    gene, tax, mog, op = p[0], p[1], p[5], p[6]
    # if op in op_mogs and tax == '511145':
    if op in operons:
        if gene in coords:
            c = coords[gene]+[gene, mog, op]
            a[tax].append(c)
g.close()

# сортируем по первой координате
for tax in a:
    a[tax].sort(key=lambda x: x[0])

'''print len(taxids)
if '511145' in taxids:
    print 'yes'
else:
    print 'no'
print 'length of a', len(a)
for i in a['511145']:
    print i'''

# извлекаем скоры по координатам
for tax in a:
    strand = {'+':'false', '-':'true'}
    # s = open('/mnt/mapr/user/chezoya/genomeScores/'+tax)
    s = open(tax)
    n = 0
    scores = []
    for line in s:
        if n < len(a[tax]):
            p = line.rstrip('\n').split()
            if strand[p[5]] == a[tax][n][2] and int(p[1]) >= a[tax][n][0] and int(p[2]) <= a[tax][n][1]:
                scores.append(float(p[4]))
            if int(p[1]) > a[tax][n][1]:
                if scores:
                    a[tax][n].append(min(scores))
                else:
                    a[tax][n].append(0)
                n += 1
                scores = []
        else:
            break
    s.close()

'''print
for i in a['511145']:
    print i'''

# теперь нет необходимости держать списки про каждый из генов, сделаем словарь
# запишем все найденные скоры в файлы по taxID (чтобы_было)
for tax in a:
    o = open('interScores/'+tax, 'a')
    o.write('locusID\tMOG\ttuID\tscore\n')
    for gene in a[tax]:
        o.write(gene[3]+'\t'+gene[4]+'\t'+gene[5]+'\t'+str(gene[6])+'\n')
    a[tax] = {x[3]:[x[4],x[6]] for x in a[tax]}

# создаём словарь предыдущих в опреоне генов
pio = {}
k = open('../info/PIO_all')
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
            o = open('operonScores/'+mog, 'a')
            o.write(tax+'\t'+str(min(s))+'\n')
            o.close()
