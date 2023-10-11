"""creates a text file in which x=mutation position, y=frequency from a mutation export .csv file
(downloaded from https://cancer.sanger.ac.uk/cosmic by searching for gene/ variants/ export:csv ) 
"""

with open (input('Enter .csv mutation export file pathname: '),'r') as fin:
    """Prompt user for csv file, open it, strip its lines for '\n' and 
    include them into a list
    """
#with open ('/Users/federicopozzani/Desktop/Lab/FLCN/clustering/6ulg_1/v2whole/sh_cluster/mut_freq/FLCN.csv') as fin:
    mylist= list (fin)
    mutlist = [line.rstrip ('\n') for line in mylist]


def selectmissense ():
    miss_list= []
    for ele in mutlist:
        if 'Missense' in ele:
            miss_list.append (ele)

    return (miss_list)

misslist= selectmissense ()

#print (misslist)


def cutandkeep():
    mutfreq=[]
    for ele in misslist:
        delimiter=','
        t=ele.split(delimiter)
        
        a=(t[2][3:-1]), t[4]
        mutfreq.append(a)
    return mutfreq
       
mutfreq= cutandkeep ()

f = open("output.txt", "w")
for row in mutfreq:
    print(" ".join(row), file=f)
f.close()


