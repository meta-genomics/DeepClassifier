import datetime
import os.path
import gzip 


kmerSize = 5

def extractKmers(p, k=kmerSize):
    kmers = [0]*(len(p)-kmerSize)
    for i in range(len(p)-kmerSize):
        kmer = 0
        for j in range(i,i+kmerSize):
            if p[j] == 'A':
                c = 0
                kmer += c * 4** (j-i)
                continue
            elif p[j] == 'C':
                c = 1
                kmer += c * 4** (j-i)
                continue
            elif p[j] == 'G':
                c = 2
                kmer += c * 4** (j-i)
                continue
            elif p[j] == 'T':
                c =3
                kmer += c * 4** (j-i)
                continue
        kmers[i] = kmer
    return kmers

dataPath = '/home/fonty/data/genomics/unam/maribel'
out = open(os.path.join(dataPath,f"kmerCounts.csv.gz"),'wt')
runs = open(os.path.join(dataPath,"fileList.txt"),'rt')

runNumber = 0

for runName in runs:
    runNumber += 1
    runName = runName.replace('\n','')
    r1Name = runName + '_R1.fastq.gz'
    r2Name = runName + '_R2.fastq.gz'

    r1 = gzip.open(os.path.join(dataPath,r1Name),'rt')
    r2 = gzip.open(os.path.join(dataPath,r2Name),'rt')

    print(f"----------------------------------")
    print(f"---- Run {runNumber} file {runName} ------")

    start = datetime.datetime.now()
    count = 0

    for p,q in zip(r1,r2):
        kmerCounts = [0] * 4**kmerSize
        count += 1
        if count % 100_000 == 0:
            currTime = datetime.datetime.now()
            timelapse = currTime-start
            progress = count / 1_200_001
            print(f"runtime:\t{timelapse}")
#            print(f"count:\t{count}")
            print(f"progress:\t{progress*100}")        
            print ( f"time to finish: \t {(timelapse / progress) - timelapse}")
            #break
        # skip read name
        # get read
        p = r1.readline()
        # comment line
        r1.readline()
        # quality line
        r1.readline()

        
        q = r2.readline()
        r2.readline()
        r2.readline()

        
        kmers = extractKmers(p,kmerSize)
        for k in kmers:
            kmerCounts[k] += 1
        kmers = extractKmers(q,kmerSize)
        for k in kmers:
            kmerCounts[k] += 1

        out.write(f"{runName}")
        for k in range(len(kmerCounts)):
            out.write(f",{kmerCounts[k]}")
        out.write("\n")


        