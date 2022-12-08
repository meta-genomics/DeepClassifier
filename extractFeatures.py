import datetime
import os.path
import gzip 

NUMBER_OF_NUCLEOTIDES = 4
NUCLEOTIDES_DICT = {'A':0,'a':0,'T':1,'t':1,'G':2,'g':2,'C':3,'c':3}
ENCODED_NUCLEOTIDES_DICT = {'A':[1,0,0,0],'a':[1,0,0,0],'T':[0,1,0,0],'t':[0,1,0,0],
                            'G':[0,0,1,0],'g':[0,0,1,0],'C':[0,0,0,1],'c':[0,0,0,1]}

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

def FastqToDataFrame(fileList, test = False, encoding='onehot', testSize=100_000):
    start = datetime.datetime.now()
    classNum = 0
    sequences = list()
    classes = list()

    for fileNameR1 in fileList:
        currTime = datetime.datetime.now()
        timelapse = currTime-start    
        print(f"runtime:\t{timelapse}\t{fileNameR1}")

        classNum += 1
        if "_R1.fastq" in fileNameR1:
            fileNameR2 = fileNameR1.replace("_R1.fastq","_R2.fastq")
        r1 = gzip.open(fileNameR1,'rt')
        r2 = gzip.open(fileNameR2,'rt')
        testCount = 0

        for p,q in zip(r1,r2):
            if test:
                testCount += 1
                if testCount > testSize:
                    break
        # skip read name
        # get read
            p = r1.readline().replace('\n','')
            q = r2.readline().replace('\n','')
        # comment line
            r1.readline()
            r2.readline()
        # quality line
            r1.readline()
            r2.readline()

            sequences.append(list(p))
            sequences.append(list(q))
            classes.append(classNum)
            classes.append(classNum)
    numClasses = classNum
    classesVector = list()
    for c in classes:
        vector = [0]*numClasses
        vector[c-1] = 1
        classesVector.append(vector)

    encodedSequences = list()
    for s in sequences:
        encodedSequence = list()
        for nt in s:
            #encoded = [0] * NUMBER_OF_NUCLEOTIDES
            if encoding == 'token':
                encodedSequence.append(NUCLEOTIDES_DICT[nt])
            elif encoding == 'vector':
                encodedSequence.append(ENCODED_NUCLEOTIDES_DICT[nt])
        # one hot sequential
            else:
                encodedSequence.extend(ENCODED_NUCLEOTIDES_DICT[nt])
        encodedSequences.append(encodedSequence)
    return (encodedSequences, classesVector, numClasses) #pd.DataFrame(list(zip( sequences, classess)), columns=[list(range(len(sequences[0]))),'class'])
            

            


        

if __name__ == "__main__":

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
