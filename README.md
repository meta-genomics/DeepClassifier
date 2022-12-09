# DeepClassifier

This code provides the proof of concept to predict taxonomic assignment of metagenomic sequences using Neural Networks.
The code uses tensor flow for the Neural Networks.

Command to run

python neuralNetwork.py --fileList fileList.txt
--fileList Obligatory parameter with the list of fastq files, paired end is assumed and only "R1" must be provided, "R2" file is assumed to exist.
--path Optional parameter for the fastq files directory, default is root folder
--epochs (30) Optional epochs parameter for TensorFlow training
--batchSize (60) Optional batchSize parameter for Tensor flow training
