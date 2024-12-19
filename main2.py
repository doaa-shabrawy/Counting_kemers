#####################################################
######## WELCOME TO THE KMER FINDER PROGRAM #########
#####################################################

##################### IMPORTS #######################
import os, sys, time
from collections import Counter
from joblib import Parallel, delayed
import heapq
from tkinter import *
from tkinter import ttk 
from tkinter.simpledialog import askstring, askinteger
from tkinter import filedialog

##################### FUNCTIONS #####################
#### FASTA INDEXING FUNCTION ###
def indexfasta(filename):
    # Open file
    infile = open(filename, 'rb')

    # Set chunksize
    chunksize = 1024 * 1024
    filepos = 0
    headstart = list()
    headend = list()

    # Read chunck size
    while True:
        content = infile.read(chunksize)
        if len(content) == 0:
            break

        # Find Header Start
        chunkpos = 0
        while chunkpos != -1:
            chunkpos = content.find(b'>', chunkpos)
            if chunkpos != -1:
                headstart.append(chunkpos + filepos)
                chunkpos += 1

        # Find Header End
        for i in range(len(headend), len(headstart)):
            chunkpos = max(0, headstart[i] - filepos)
            chunkpos = content.find(b'\n', chunkpos)
            if chunkpos != -1:
                headend.append(chunkpos + filepos)
        filepos += len(content)
    infile.close()

    # Eliminating wrong headers due to extra > in header line
    for i in range(len(headstart) - 1, 0, -1):
        if headend[i] == headend[i - 1]:
            del headstart[i]
            del headend[i]
    headstart.append(filepos)
    fastaindex = list()
    for i in range(len(headend)):
        seq_start = headend[i] + 1
        seq_end = headstart[i + 1] - 1
        fastaindex.append((seq_start, seq_end, seq_end - seq_start))

    # Load Balancing
    # fastaindex = sorted(fastaindex, key=lambda x: x[2], reverse=False)

    return fastaindex


#### ENTRY INDEXING FUNCTION ####
def indexsequence(seq):
    pointer = 0
    seqindex = list()

    while len(seq) > pointer:

        # Find start of seq
        potenstart = [seq.find(b'a', pointer), seq.find(b't', pointer), seq.find(b'c', pointer),
                      seq.find(b'g', pointer)]

        realstart = min(potenstart)
        if realstart == -1:
            # happens rarely, so slow code is ok
            potenstart = [i for i in potenstart if i > -1]
            if len(potenstart) == 0:
                break
            realstart = min(potenstart)
        realend = seq.find(b'N', realstart)
        if realend == -1:
            realend = len(seq)
        seqindex.append((realstart, realend))
        pointer = realend
    return seqindex


#### KMER FINDING FUNCITON ####
def find_kmers(fasta, idx):
    # Read sequence
    infile = open(fasta, 'rb')
    infile.seek(idx[0])
    seq = infile.read(idx[1] - idx[0] + 1).translate(transtable, b'\r\n\t ')
    infile.close()
    subdict = dict()

    # Index sequence
    seqindex = indexsequence(seq)

    # Slide through sequences and add to dictionary removing unwanted bases
    for start, stop in seqindex:
        for i in range(start, stop - kmer_len + 1):
            kmer = seq[i:i + kmer_len]
            if kmer not in subdict:
                subdict[kmer] = 1
            else:
                subdict[kmer] += 1

    return subdict


if __name__ == '__main__':

    # SET VARIABLES
    window=Tk()
    file=filedialog.askopenfilename(parent=window)
    kmer_len= askinteger('K-Mer Counting Program' ,'Please enter the lengh of kemer',parent=window)
    
    transtable = bytes.maketrans(b'ATCGMRYKVHDBWmrykvhdbxnsw', b'atcgNNNNNNNNNNNNNNNNNNNNN')
    n_worker = os.cpu_count()

    # INDEXING
    start = time.time()
    my_indexes = indexfasta(file)
    index_time = time.time() - start

    # FIND KMER
    start = time.time()
    results = Parallel(n_jobs=5)(delayed(find_kmers)(file, z) for z in my_indexes)
    search_time = time.time() - start

    # MERGE DICTS
    final_dict = dict()
    start = time.time()
    for r in results:
        for merkey in r.keys():
            if merkey in final_dict:
                final_dict[merkey] += r[merkey]
            else:
                final_dict[merkey] = r[merkey]
    dict_sum_time = time.time() - start
    del results
    del r

    # displaying sorted list of kemers

    targetfile=open(file)
    readfile=targetfile.read()
    tseq ="".join(readfile.split())
    targetfile.close()
    def cmers (tseq,kmer_len):
        kfreq={}
        for i in range(0,len(tseq)-kmer_len+1):
            kmer=tseq[i:i+kmer_len]
            if kmer in kfreq:
                kfreq[kmer] +=1
            else:
                kfreq[kmer]=1
        return kfreq

    rf=cmers(tseq,kmer_len)
    # sorted alphabitical
    listrf_alpha=[[seq,rf[seq]]for seq in rf]
    listrf_alpha.sort()
    #sorted numircal
    textfile =open("a_file.fsa", "w")
    for element in listrf_alpha:
            textfile.write(str(element) + ',')
    textfile.close()
       
        
    def display5():
        listrf=[[rf[seq],seq]for seq in rf]
        listrf.sort(reverse=True)
        maxn=askinteger('K-Mer Counting Program' ,'enter the number of top n of kemers')
        top_n=heapq.nlargest(maxn, listrf)
        j=240
        for i in range(0,len(top_n)):
            label=Label(window,text=top_n[i],font=('Great Vibes',13),bg="#33A076")
            label.place(x=10,y=j)
            j+=25
        
 
    def display1():
        
        e1=round(index_time, 3)
        label1.config(text=e1)
        
    def display2():
        
        e1=round(search_time, 3)
        label2.config(text=e1)
        
    def display3():
        
        e1=round(dict_sum_time, 3)
        label3.config(text=e1)
        
    def display4():
        
        e1= round(index_time + search_time + dict_sum_time, 3)
        label4.config(text=e1)
        
   
        
        
    window.title("k-mer counting program")
    #creating Labels
    label=Label(window,text='K-Mer Counting Program',font=('Courier',30,'bold'),bg="#3e3e42",fg="white")
    label.pack(side=TOP,fill=X)
    label=Label(window,text='',font=('arial',25,'bold'),bg="#3e3e42")
    label.pack(side=BOTTOM,fill=X)
    label=Label(window,text='RESULTS FOR KMER FINDER:',font=('Great Vibes',13) ,bg="#33A076")
    label.place(x=30,y=60)    
       
    label1=Label(window,text='',font=('Great Vibes',13),bg="#33A076")
    label1.place(x=40,y=130) 
    
    label2=Label(window,text="",font=('Great Vibes',13),bg="#33A076")
    label2.place(x=240,y=130) 
    
    label3=Label(window,text="",font=('Great Vibes',13),bg="#33A076")
    label3.place(x=440,y=130) 
    
    label4=Label(window,text="",font=('Great Vibes',13),bg="#33A076")
    label4.place(x=340,y=210)
    
    label5=Label(window,text="",font=('Great Vibes',13),bg="#33A076")
    label5.place(x=0,y=240)
    
    label5=Label(window,text="*NOTE :"+"\n"+"All kmers and there occurences are in"+"\n"+" a file named 'a_file.fsa' presistant in "+"\n"+"the same directory of the program.",font=('Great Vibes',13),bg="#33A076")
    label5.place(x=310,y=240)
    
    
    btn=ttk.Button(window,text='Indexing time:',command=display1)
    btn.place(x=10,y=100,width=125,height=30)
    
    btn=ttk.Button(window,text="Finding kmers time:",command=display2)
    btn.place(x=210,y=100,width=125,height=30)
    
    btn=ttk.Button(window,text="Merging dicts time:",command=display3 )
    btn.place(x=410,y=100,width=125,height=30)
    
    btn=ttk.Button(window,text="TOTAL:",command=display4)
    btn.place(x=310,y=180,width=125,height=30)
    
    btn=ttk.Button(window,text='TOP_N of k-mer occurences:',command=display5)
    btn.place(x=10,y=180,width=170,height=30)
    
    
    window.geometry('600x400')
    window.configure(bg="#33A076")
    window.resizable(True,True)
    window.mainloop()
