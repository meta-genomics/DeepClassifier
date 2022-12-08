import tensorflow as tf
from tensorflow import keras
from keras import layers
from extractFeatures import *
import sklearn.metrics as metrics
import numpy as np
from keras.utils import np_utils
import random

epochs = 20
batchSize = 64

path = '/mnt/d/data/genomics/unam/maribel'
fileList = open( os.path.join(path,"fileList.txt"))
fileNames = list()
for line in fileList:
    fileNames.append(os.path.join(path,line.replace('\n','')))

#sequences, classes, numClasses = FastqToDataFrame(fileNames, test=True, encoding='token', testSize=1_00)
sequences, classes, numClasses = FastqToDataFrame(fileNames, test=True, encoding='onehot', testSize=1_00)
lengthRead = len(sequences[0])

#suffle samples
temp = list(zip(sequences, classes))
random.shuffle(temp)
sequences, classes = zip(*temp)
sequences, classes = list(sequences), list(classes)


N = len(sequences)
N90 = int(.9*N)

x_train = sequences[0: N90]
x_test = sequences[N90:-1]

y_train = classes[0:N90]
y_test = classes[N90:-1]

model = keras.models.Sequential()
model.add(layers.InputLayer(input_shape=(lengthRead,1), name='Input_Layer'))
model.add(layers.LSTM(128, return_sequences=False, activation='relu'))
#model.add(layers.SimpleRNN(128, return_sequences=False, activation='relu'))
model.add(layers.Dense(numClasses, activation='Softmax', name='Output_Layer'))
model.summary()
model.compile(optimizer='adam',
            loss='categorical_crossentropy',#sparse_categorical_crossentropy    categorical_crossentropy
            metrics=['accuracy'])
model.fit(x = x_train, 
        y = y_train,##Cuando haya otra base, se modifica ese 200
        batch_size=batchSize,
        epochs=epochs,
        verbose=2,shuffle=True)


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
