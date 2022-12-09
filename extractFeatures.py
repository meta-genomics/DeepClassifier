import datetime
import os.path
import gzip 

NUMBER_OF_NUCLEOTIDES = 4
NUCLEOTIDES_DICT = {'A':0,'a':0,'C':1,'c':1,'G':2,'g':2,'T':3,'t':3}
ENCODED_NUCLEOTIDES_DICT = {'A':[1,0,0,0],'a':[1,0,0,0],'C':[0,1,0,0],'c':[0,1,0,0],
                            'G':[0,0,1,0],'g':[0,0,1,0],'T':[0,0,0,1],'t':[0,0,0,1]}

kmerSize = 5
# get kmers from a sequence of nucleotides codified in base 4
# p is the sequence and k the size of the kmers
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

# Three types of encoding are provided "onehot" is serial one hot, "vector" is one hot as a vector and "token" are nucleotides as integers
# File list has all the files with only the R1 paired end and will look for the same file with "R2"
# If test is True, the maximum number of reads processed is given by testSize
def FastqFilesToFeatures(fileList, test = False, encoding='onehot', testSize=100_000):
    start = datetime.datetime.now()
    classNum = 0
    sequences = list()
    classes = list()

    for fileNameR1 in fileList:
        currTime = datetime.datetime.now()
        timelapse = currTime-start    
        print(f"runtime:\t{timelapse}\t{fileNameR1}")
        classNum += 1
    # Read the R2 paired end file
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
        # skip comment line
            r1.readline()
            r2.readline()
        # skip quality line
            r1.readline()
            r2.readline()
            sequences.append(list(p))
            sequences.append(list(q))
            classes.append(classNum)
            classes.append(classNum)
    numClasses = classNum
    classesVector = list()
# one hot encoding de la clase 
    for c in classes:
        vector = [0]*numClasses
        vector[c-1] = 1
        classesVector.append(vector)

# encoding de la secuencia
    encodedSequences = list()
    for s in sequences:
        encodedSequence = list()
        for nt in s:
        # token A = 0, T = 1, G = 2, C =3 
            if encoding == 'token':
                encodedSequence.append(NUCLEOTIDES_DICT[nt])
        # vector A =[ 1,0,0,0], T = [0,1,0,0], ...                
            elif encoding == 'vector':
                encodedSequence.append(ENCODED_NUCLEOTIDES_DICT[nt])
        # one hot sequential AT = [1,0,0,0,0,1,0,0]
            else:
                encodedSequence.extend(ENCODED_NUCLEOTIDES_DICT[nt])
        encodedSequences.append(encodedSequence)
    return (encodedSequences, classesVector, numClasses) #pd.DataFrame(list(zip( sequences, classess)), columns=[list(range(len(sequences[0]))),'class'])
    