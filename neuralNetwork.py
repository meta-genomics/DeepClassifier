import tensorflow as tf
from tensorflow import keras
from keras import layers
from extractFeatures import *
import sklearn.metrics as metrics
import random
import os
import argparse

parser = argparse.ArgumentParser(description='Taxonomic assignation using neural networks.')
parser.add_argument('--epochs', metavar='-e', type=int,
                    help='Number of epochs for tensor flow training')

parser.add_argument('--batchSize', metavar='-b', type=int, 
                    help='Batch size for tensor flow training')

parser.add_argument('--path', metavar='-p',
                    help='Path for read filest')

parser.add_argument('--fileList', metavar='-f',
                    help='Read file list including only R1 paired end')


args = parser.parse_args()

epochs = args.epochs
batchSize = args.batchSize
path = args.path
fileList = args.fileList

if epochs == None:
    epochs = 30

if batchSize == None:
    batchSize = 60

if path == None:
    path = ''

if fileList == None:
    print("fileList is an obligatory parameter")
    exit()

print(f"epochs:   \t{epochs}")
print(f"batchSize:\t{batchSize}")
print(f"path:     \t{path}")
print(f"fileList: \t{fileList}")
#print(args)
exit()
          

# parameters for tensor flow
epochs = 30
batchSize = 64
# path where reads are stored
path = '/mnt/d/data/genomics/unam/maribel'
readFileList = 'fileList.txt'

# list of files to be processed
fileList = open( os.path.join(path,readFileList))
fileNames = list()
# pre pend the path to all the files
for line in fileList:
    fileNames.append(os.path.join(path,line.replace('\n','')))

# load all the sequences and their classes 
sequences, classes, numClasses = FastqFilesToFeatures(fileNames, test=True, encoding='onehot', testSize=10_000)
lengthRead = len(sequences[0])

# suffle samples 
temp = list(zip(sequences, classes))
random.shuffle(temp)
sequences, classes = zip(*temp)
sequences, classes = list(sequences), list(classes)

# split train test data
N = len(sequences)
N90 = int(.9*N)

x_train = sequences[0: N90]
x_test = sequences[N90:-1]

y_train = classes[0:N90]
y_test = classes[N90:-1]

# generate tensor flow model
model = keras.models.Sequential()
model.add(layers.InputLayer(input_shape=(lengthRead,1), name='Input_Layer'))
model.add(layers.LSTM(512, return_sequences=False, activation='relu'))
#model.add(layers.SimpleRNN(128, return_sequences=False, activation='relu'))
model.add(layers.Dense(numClasses, activation='Softmax', name='Output_Layer'))
model.summary()
model.compile(optimizer='adam',
            loss='categorical_crossentropy',#sparse_categorical_crossentropy    categorical_crossentropy
            metrics=['accuracy'])
model.fit(x = x_train, 
        y = y_train,
        batch_size=batchSize,
        epochs=epochs,
        verbose=2,shuffle=True)



# get performance metrics 
#model.evaluate(x_test, y_test, batch_size=batchSize, verbose=2)
test_predicted_labels_raw = model.predict(x_test)
test_true_labels      = np.argmax(y_test,axis=1)
test_predicted_labels = np.argmax(test_predicted_labels_raw,axis=1)
accuracy = metrics.accuracy_score(test_true_labels,test_predicted_labels)
precision = metrics.precision_score(test_true_labels,test_predicted_labels,average='weighted')
recall = metrics.recall_score(test_true_labels,test_predicted_labels,average='weighted')
f1 = metrics.f1_score(test_true_labels,test_predicted_labels,average='weighted')


print(accuracy)
print(precision)
print(recall)
print(f1)
