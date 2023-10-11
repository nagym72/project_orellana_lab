#!/usr/bin/env python
# coding: utf-8

# Retrieves PDB IDs from Uniprot Accession Number
# Interactively prompts for input:
#   - Uniprot ID
#   - Range of perc identity: 100/90/50
# Outputs files in the working directory

import requests, sys, json

# Documentation: https://www.uniprot.org/help/api
WEBSITE_API = "https://rest.uniprot.org/"

# Documentation: https://www.ebi.ac.uk/proteins/api/doc/
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

# Helper function to download data
def get_url(url, **kwargs):
  response = requests.get(url, **kwargs);

  if not response.ok:
    print(response.text)
    response.raise_for_status()
    sys.exit()

  return response


# In[13]:


accession = input ('Enter Uniprot accession number:\n')
idperc = input ('Enter identity % (100/90/50):\n')
#accession = 'P61626'
#idperc = '50'


# In[14]:


def createfiles(acc):
    if idperc == '100':
        c = 1
        r = get_url(f"{WEBSITE_API}uniprotkb/stream?format=txt&query=%28uniref_cluster_100%3AUniRef100_{acc}%29")
        out= (r.text)
        fp = open('ide_100%.txt', 'w')
        fp.write(out)
        fp.close()
        return (c)
    
    if idperc == '90':
        c = 2
        o = get_url(f"{WEBSITE_API}uniprotkb/stream?format=txt&query=%28uniref_cluster_90%3AUniRef90_{acc}%29")
        whole= (o.text)
        rp = open('ide_90%.txt', 'w')
        rp.write(whole)
        rp.close()
        r = get_url(f"{WEBSITE_API}uniprotkb/stream?format=txt&query=%28uniref_cluster_100%3AUniRef100_{acc}%29")
        top100= (r.text)
        op = open('ide_100%.txt', 'w')
        op.write(top100)
        op.close()
        return (c)
        
    if idperc == '50':
        c = 3
        p = get_url(f"{WEBSITE_API}uniprotkb/stream?format=txt&query=%28uniref_cluster_50%3AUniRef50_{acc}%29")
        whole= (p.text)
        pp = open('ide_50%.txt', 'w')
        pp.write(whole)
        pp.close()
        o = get_url(f"{WEBSITE_API}uniprotkb/stream?format=txt&query=%28uniref_cluster_90%3AUniRef90_{acc}%29")
        top90= (o.text)
        rp = open('ide_90%.txt', 'w')
        rp.write(top90)
        rp.close()
        r = get_url(f"{WEBSITE_API}uniprotkb/stream?format=txt&query=%28uniref_cluster_100%3AUniRef100_{acc}%29")
        top100= (r.text)
        op = open('ide_100%.txt', 'w')
        op.write(top100)
        op.close()
        return (c)
    
c = createfiles(accession)

print ( str(c+1)+' files created')


# In[15]:


def readfile(filein):
    #with open ('./output.txt') as fin:
    with open (filein) as fin:

        mylist= list (fin)
        outp = [line.rstrip ('\n') for line in mylist]
    return outp

#print (readfile('./ide_100%.txt'))


# In[33]:


def idpdb(fin):
    l=[]
    for line in fin:
        if line[:2] == 'ID':
            prot=(line[5:])
            #print (prot)


        if line[:9] == 'DR   PDB;':
            pdb=(line[10:])
            #print (pdb)
            #id2pdb[prot]=pdb
            l.append(pdb+ (' '*(58 - len(line))) + prot)
            #if len(line) > 40:
              #l.append(pdb+'\t'+prot)
            #if len(line) < 36:
              #l.append(pdb+'\t\t\t'+prot)
            #else:
              #l.append(pdb+'\t\t'+prot)

    return(l)

#for line in (idpdb(readfile('./ide_50%.txt'))):
    #print (line)

#fp = open('prova.txt', 'w')
#fp.write(idpdb(readfile('./ide_100%.txt')))
#fp.close()    


# In[34]:


def sort():
    l=[]
    if idperc == '100':
        for line in idpdb(readfile('./ide_100%.txt')):
            l.append(line+'\t 100%')
            
    if idperc == '90':
        for line in idpdb(readfile('./ide_100%.txt')):
            l.append(line+'\t 100%')
        for line in idpdb(readfile('./ide_90%.txt')):
            if line not in idpdb(readfile('./ide_100%.txt')):
                l.append(line+'\t 90-100%')
    if idperc == '50':
        for line in idpdb(readfile('./ide_100%.txt')):
            l.append(line+'\t 100%')
        for line in idpdb(readfile('./ide_90%.txt')):
            if line not in idpdb(readfile('./ide_100%.txt')):
                l.append(line+'\t 90-100%')
        for line in idpdb(readfile('./ide_50%.txt')):
            if line not in idpdb(readfile('./ide_100%.txt')) and line not in idpdb(readfile('./ide_90%.txt')):
                l.append(line+'\t 50-100%')
            
    return (l)
    

#for line in sort():
    #print (line)
f = open("pdbs.txt", "w")
for row in sort():
    print (row, file=f)
f.close()
print ('output file created')
print ('job done!')


# In[ ]:




