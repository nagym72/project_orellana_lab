""" Takes a consurf summary .txt file as imput and cleans it into a scores.dat file
in which x = residue, y = consurf score
"""

with open (input('Enter Consurf Summary file pathname: '),'r') as fin:
    #"""Prompt user for consurf summary .txt file, open it, strip its lines for '\n' and 
    #include them into a list
    #"""

#with open ('/Users/federicopozzani/Desktop/Lab/FLCN/clustering/6ulg_1/v2whole/sh_cluster/clean_consurf/6ULGL_consurf_summary.txt') as fin:
    

    mylist= list (fin)
    mutlist = [line.rstrip ('\n') for line in mylist]
    
del mutlist[:15]
del mutlist[-4:]


def delim():
    split=[]
    
    for ele in mutlist:
        delimiter='\t'
        t=ele.split(delimiter)
        if '-' not in t[2]:
            a=t[0].strip(' '), t[5].strip(' *')
            split.append(a)
    return split


f = open("scores.dat", "w")
for row in delim():
    print (" ".join(row), file=f)
    
f.close()
    
