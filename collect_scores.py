# -*- coding: utf-8 -*-
# скрипт берёт: файл-список интересующих MOG; файл про все гены, файл с координатами апстримов, файлы со скорами
# скрипт выдаёт для каждого MOG: файл со скором в кажом геноме
import os.path
import csv

# делаем словарь с интересующими нас MOG и их именами
mog_names = {}
f = open('rgenes/rprots.txt')
fr = csv.DictReader(f, delimiter='\t')
for line in fr:
    mog_names[line['MOG']] = line['protein']
f.close()
print 'mog_names', len(mog_names)

# делаем словарь с координатами апстримов для всех генов
coords = {}
h = open('info/upstreams-350+50-80')
hr = csv.DictReader(h, delimiter='\t')
for line in hr:
    coords[line['MO']] = [int(line['upsBegin']), int(line['upsEnd']), line['thisGeneIsCompl'], int(line['upsLength'])]
h.close()
print 'coords', len(coords)

# выписываем номера оперонов, в которых лежат нужные MOG, и номера таксонов
taxids = []
operons = []
g = open('info/allBacteriaGenesInfo2.csv')
gr = csv.DictReader(g, delimiter='\t')
for line in gr:
    # if line['mogId'] in mog_names and line['taxonomyId'] == '511145':
    if line['mogId'] in mog_names:
        if line['tuId'] not in operons:
            operons.append(line['tuId'])
        if line['taxonomyId'] not in taxids:
            taxids.append(line['taxonomyId'])
g.close()
print 'taxids', len(taxids)
print 'operons', len(operons)

# создаём словарь геном:список_координат всех генов из оперонов, попавших в operons
# такие списки можно будет отсортировать, чтобы искать скоры
a = {x:[] for x in taxids}
q = open('info/allBacteriaGenesInfo2.csv')
qr = csv.DictReader(q, delimiter='\t')
for line in qr:
    if line['tuId'] in operons:
        if line['\xef\xbb\xbflocusId'] in coords:
            c = coords[line['\xef\xbb\xbflocusId']]+[line['\xef\xbb\xbflocusId'], line['mogId'], line['tuId']]
            a[line['taxonomyId']].append(c)
q.close()
print 'a',len(a)

# сортируем по первой координате
for tax in a:
    a[tax].sort(key=lambda x: x[0])

# извлекаем скоры по координатам
direction = {'+':'false', '-':'true'}
unavailable = []
for tax in a:
    fname = '/mnt/mapr/user/chezoya/genomeScores/'+tax
    # fname = 'rgenes/'+tax
    if os.path.isfile(fname):
        s = open(fname)
        sr = csv.DictReader(s, fieldnames=['chr','start','end','energy','score','strand'], delimiter='\t')
        n = 0
        scores = []
        for line in sr:
            if n < len(a[tax]):
                gene_start, gene_end, gene_strand = a[tax][n][0], a[tax][n][1], a[tax][n][2]
                if direction[line['strand']] == gene_strand and gene_start <= int(line['start'])<= int(line['end']) <= gene_end:
                    scores.append(float(line['score']))
                if int(line['end']) > gene_end:
                    a[tax][n].append(min(scores+[0]))
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

print a['511145']
# теперь нет необходимости держать списки про каждый из генов, сделаем словарь
# запишем все найденные скоры в файлы по taxID (чтобы_было)
for tax in a:
    o = open('rgenes/interScores/'+tax, 'a')
    o.write('locusID\tMOG\ttuID\tupstLen\tscore\n')
    for gene in a[tax]:
        if len(gene)!=8:
            print 'tax', tax
            print 'gene', gene
            print a[tax]
        o.write('\t'.join(gene[4:7])+'\t'+str(gene[3])+'\t'+str(gene[7])+'\n')
    a[tax] = {x[4]:[x[5],x[7],x[3]] for x in a[tax]}
    o.close()

# создаём словарь предыдущих в опреоне генов
pio = {}
k = open('info/PIO_all')
kr = csv.DictReader(k, delimiter='\t')
for line in kr:
    pio[line['geneId']] = pio.get(line['geneId'], []) + [line['geneIdPiO']]
k.close()

# выбираем из скоров в опероне только те, что лежат выше гена, минимальный записываем в файл
for tax in a:
    for gene in a[tax]:
        mog = a[tax][gene][0]
        if mog in mog_names:
            s = []
            for i in pio[gene]:
                if i in a[tax]:
                    s.append(a[tax][i][1:])
            s.sort(key=lambda x: x[0])
            o = open('rgenes/operonScores/'+mog, 'a')
            o.write(tax+'\t'+str(s[0][0])+'\t'+str(s[0][1])+'\n')
            o.close()
