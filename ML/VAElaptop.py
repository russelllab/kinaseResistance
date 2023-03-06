#!/usr/bin/env python
# coding: utf-8

## Develop VAE for kinase resistance predtictions
import numpy as np
import scipy as sp
import os, sys, gzip
import seaborn as sns
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn.preprocessing import MinMaxScaler
tf.random.set_seed(69)
from cls import Kinase, Mutation
import fetchData

PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

exceptions= ['Q9Y2K2', 'Q15303', 'Q9UIK4', 'P33981', 'P35916',
             'Q99683', 'Q8IVW4', 'Q9H792', 'Q9P286', 'Q86Z02',
             'Q8TF76', 'Q96L34', 'Q13308', 'Q9UK32', 'Q15772',
             'P51617', 'Q9Y3S1', 'Q9C098', 'Q6VAB6', 'P21127',
             'Q13557', 'Q6ZMQ8', 'Q6P0Q8', 'Q8IZE3', 'P51957',
             'O60229', 'Q96RG2', 'Q5VST9', 'Q8WZ42', 'O75962',
             'O95835', 'Q13535']


kinases = {}

fetchData.fetchFasta(kinases, Kinase)
fetchData.fetchGroup(kinases, Kinase)
# print (kinases['P00533'].group)
hmmPkinase = fetchData.fetchPkinaseHMM() # hmmPosition > AA > bit-score
# print (hmmPkinase[30]['K'])
fetchData.fetchHmmsearch(kinases, Kinase)
# fetchData.dsspScores(kinases, Kinase)
# # print (kinases['Q9NYV4'].burr[3])
# # print (kinases['Q92772'].dihedral)
# fetchData.iupredScores(kinases, Kinase)
# fetchData.homologyScores(kinases, Kinase)

# #print (kinases['Q9NYV4'].mechismo)
# data = []
# for acc in kinases:
#     #print (kinases[acc].domains)
#     for domainNum in kinases[acc].domains:
#         #print (domainNum)
#         data.append(len(kinases[acc].domains[domainNum]))

# #print (data)
# df = pd.DataFrame(data = data, columns=['Length'])
# sns.histplot(data=df, x="Length")

def fetchPkinase(acc, domainNum):
    '''
    A function to take acc and alignment cutoff values to return
    the HMM bitscores
    '''
    row = []
    for hmmPosition in range(1,265):
        if hmmPosition in kinases[acc].domains[domainNum]:
            SeqPosition = kinases[acc].domains[domainNum][hmmPosition]
            aa = kinases[acc].fasta[SeqPosition-1]
            value = float(hmmPkinase[hmmPosition][aa])
        else:
            value = 3

        row.append(value)
    #print (len(row))
    kinases[acc].hmm[domainNum] = row
    return row

def fetchStrucFeat(acc, domainNum):
    data = []
    df = pd.DataFrame()
    for dic, name in zip([
                        kinases[acc].dihedral,
                        kinases[acc].sec,
                        kinases[acc].burr,
                        kinases[acc].access,
                        kinases[acc].iupred,
                        kinases[acc].mechismo
                        ],[
                        'dihedral',
                        'access',
                        'burr',
                        'sec',
                        'iupred',
                        'mechismo'
                        ]):
        row = []
        print (kinases[acc].dihedral)
        sys.exit()
        for hmmPosition in range(1,265):
            if hmmPosition in kinases[acc].domains[domainNum]:
                SeqPosition = kinases[acc].domains[domainNum][hmmPosition]
                residue = kinases[acc].fasta[SeqPosition-1]
                #print (SeqPosition)
                try:
                    value = dic[SeqPosition][residue]
                except:
                    print (acc, SeqPosition, len(dic), residue)
                    sys.exit()
            else:
                value = 3
            row.append(value)
        #f = pd.DataFrame(data, columns=AA)
        df[name] = row
    #print (df)

    return df

def fetchSeqFeat(acc, domainNum):
    data = []
    df = pd.DataFrame()
    for dic, name in zip([
                        kinases[acc].orthologs,
                        kinases[acc].exclParalogs,
                        kinases[acc].specParalogs,
                        kinases[acc].allHomologs,
                        kinases[acc].bpso,
                        kinases[acc].bpsh
                        ],[
                        'orthologs',
                        'exclParalogs',
                        'specParalogs',
                        'allHomologs',
                        'bpso',
                        'bpsh'
                        ]):
        row = []
        for hmmPosition in range(1,265):
            if hmmPosition in kinases[acc].domains[domainNum]:
                SeqPosition = kinases[acc].domains[domainNum][hmmPosition]
                residue = kinases[acc].fasta[SeqPosition-1]
                #print (SeqPosition)
                try:
                    #print (dic[str(SeqPosition)])
                    for mutation in dic[str(SeqPosition)]:
                        if mutation[-1] == residue:
                            value = dic[str(SeqPosition)][mutation]
                            break
                except:
                    print (acc, SeqPosition, len(dic), residue)
                    sys.exit()
            else:
                value = 3
            row.append(value)
        #f = pd.DataFrame(data, columns=AA)
        df[name] = row
    #print (df)
    return df


def oneHotEncoding(acc, domainNum):
    '''
    A function to take acc and alignment cutoff values to return
    the feature matrix (one hot encoded)
    '''
    #print (len(AA))
    data = []
    numZeros = 0
    for i in range(1,265):
        if i in kinases[acc].domains[domainNum]:
            position = kinases[acc].domains[domainNum][i]
            residue = kinases[acc].fasta[position-1]
        else:
            residue = '-'

        row = []
        for aa in AA:
            if residue == aa:
                row.append(1)
            else:
                row.append(0)
        data.append(row)
        '''
        if acc == 'Q96NX5' and i == 173:
            print (residue)
            print (row)
        '''
    data = np.array(data)
    kinases[acc].oneHotEncoding[domainNum] = data
    '''
    if data.shape == (264,len(AA)):
        #print (data.shape)
        trainData.append(data)
    '''

    #trainData = np.array(trainData)
    #print (np.stack(trainData, axis=0).shape)
    return (data, AA)

'''Map sequence to Pfam for all canonical kinases'''
seq2pfam = {}
for line in gzip.open('../data/humanKinasesHmmsearchMappings2.tsv.gz', 'rt'):
    if line[0] =='#':
        continue
    acc = line.split('\t')[0].split('|')[1]
    seqPos = str(line.split('\t')[2])
    pfamPos = str(line.replace('\n', '').split('\t')[4])
    if acc not in seq2pfam: seq2pfam[acc] = {}
    seq2pfam[acc][seqPos] = pfamPos

pkinase_act_deact_res = {'A': [], 'D': [], 'R': []}
'''Fetch act/deact mutation data'''
for line in open('../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
    if line.split()[0] == 'uniprot_name': continue
    gene = line.split('\t')[0]
    acc = line.split('\t')[1]
    wtAA = line.split('\t')[2]
    mutAA = line.split('\t')[4]
    if len(wtAA) > 1 or len(mutAA) > 1: continue
    position = str(line.split('\t')[3])
    mut_type = line.split('\t')[5]
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    mutation = wtAA + position + mutAA
    kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc)
    kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][position]
    pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)

pkinase_resistant = []
'''Fetch resistant mutation data'''
for line in gzip.open('../KA/resistant_mutations_Mar_2023.tsv.gz', 'rt'):
    if line[0] == '#': continue
    actual_gene = line.split('\t')[0]
    if '_' in actual_gene: gene = actual_gene.split('_')[0]
    else: gene = actual_gene
    acc = line.split('\t')[2]
    cosmic_mutation = line.split('\t')[1]
    wtAA = cosmic_mutation[0]
    mutAA = cosmic_mutation[-1]
    if mutAA == 'X': continue
    if len(wtAA) > 1 or len(mutAA) > 1: continue
    uniprot_position = line.split('\t')[5].replace('\n', '')
    mutation = wtAA + uniprot_position + mutAA
    mut_type = 'R'
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    if mutation not in kinases[acc].mutations:
        kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc)
        try:
            kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][uniprot_position]
            pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)
        except:
            print (acc, actual_gene, uniprot_position, mutation, cosmic_mutation)
            sys.exit()
    else: kinases[acc].mutations[mutation].mut_types.append(mut_type)
    # kinases[acc].mutations[wtAA+position+mutAA] = Mutation(wtAA+position+mutAA, mut_type)

'''Fetch PTM data'''
hmmPTM = {}
for line in open('../data/Kinase_psites4.tsv', 'r'):
    if line[0] == '#': continue
    if line.split()[2] != 'Pkinase': continue
    acc = line.split('\t')[0]
    if acc not in kinases: continue
    position = int((line.split('\t')[3].split('-')[0])[1:])
    ptm_type = line.split('\t')[3].split('-')[1]
    hmm_position = int(line.split('\t')[4])
    if ptm_type not in kinases[acc].ptm:
        kinases[acc].ptm[ptm_type] = []
    kinases[acc].ptm[ptm_type].append(position)
    if hmm_position not in hmmPTM:
        hmmPTM[hmm_position] = []
    hmmPTM[hmm_position].append(ptm_type)
    # print (acc, kinases[acc].gene, position, ptm_type)
# print (hmmPTM[141])
# sys.exit()

'''Make training matrix'''
trainMat = 'Acc\tMutation\t'
trainMat += 'hmmScoreWT\thmmScoreMUT\t'
trainMat += '\t'.join(PTM_TYPES) + '\t'
trainMat += '_pfam\t'.join(PTM_TYPES) + '_pfam\t'
trainMat += '_WT\t'.join(AA) + '_WT\t'
trainMat += '_MUT\t'.join(AA) + '_MUT\t'
# trainMat += '\t'.join(['allHomologs','exclParalogs','specParalogs','orthologs','bpso','bpsh']) + '\t'
trainMat += '_known\t'.join(['A', 'D', 'R']) + '_known\t'
trainMat += 'MUT_TYPE\n'
# print (trainMat)
# print ('_WT\t'.join(AA) + '\t')
# print (AA)
# sys.exit()
data = []
mut_types_colors = []
for acc in kinases:
    if len(kinases[acc].mutations) == 0:
        continue
    for mutation in kinases[acc].mutations:
        row = []
        mutation_obj = kinases[acc].mutations[mutation]
        position = mutation_obj.position
        mutAA = mutation_obj.mutAA
        wtAA = mutation_obj.wtAA
        mut_types = list(set(mutation_obj.mut_types))
        hmmPos, hmmScoreWT, hmmScoreMUT = fetchData.getHmmPkinaseScore(acc, wtAA, position, mutAA, kinases, hmmPkinase)
        ptm_row = fetchData.getPTMscore(acc, position, kinases, hmmPTM)
        aa_row = fetchData.getAAvector(wtAA, mutAA)
        # homology_row = fetchData.getHomologyScores(acc, wtAA, position, mutAA, kinases)
        print (
            acc +'\t'+ mutation +'\t'+ str(hmmPos) +'\t'+
            str(hmmScoreWT)+'\t' +str(hmmScoreMUT)+'\t'+ ','.join(ptm_row) + '\t' +
            ','.join(aa_row) + '\t' + '\t'.join(mut_types)
            )
        row.append(float(hmmScoreWT))
        row.append(float(hmmScoreMUT))
        row += [int(item) for item in ptm_row]
        row += [int(item) for item in aa_row]
        # row += homology_row
        adr_row = []
        for mut_type in ['A', 'D', 'R']:
            if str(hmmPos) in pkinase_act_deact_res[mut_type]: adr_row.append(1)
            else: adr_row.append(0)
        row += [int(item) for item in adr_row]
        for mut_type in mut_types:
            if mut_type == 'A':
                mut_types_colors.append('green')
            elif mut_type == 'D':
                mut_types_colors.append('red')
            else:
                mut_types_colors.append('violet')
            data.append(row)
            trainMat += acc + '\t' + mutation + '\t'
            trainMat += str(hmmScoreWT) + '\t' + str(hmmScoreMUT) + '\t'
            trainMat += '\t'.join([str(item) for item in ptm_row]) + '\t'
            trainMat += '\t'.join([str(item) for item in aa_row]) + '\t'
            # trainMat += '\t'.join([str(item) for item in homology_row]) + '\t'
            trainMat += '\t'.join([str(item) for item in adr_row]) + '\t'
            trainMat += mut_type + '\n'

gzip.open('trainData.tsv.gz', 'wt').write(trainMat)
data = np.array(data)
# scaler = MinMaxScaler()
# scaler.fit(data)
# data = scaler.transform(data)
# print (trainMat)
# sys.exit()

pca = decomposition.PCA(n_components=2)
pca.fit(data)
X = pca.transform(data)

fig = plt.figure(1, figsize=(4, 3))
plt.clf()

ax = fig.add_subplot(111)
# ax = fig.add_subplot(111, projection="3d", elev=48, azim=134)
ax.set_position([0.1, 0.1, 0.8, 0.8])


plt.cla()

# ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=mut_type_colors, cmap=plt.cm.nipy_spectral, edgecolor="k")
ax.scatter(X[:, 0], X[:, 1], c=mut_types_colors, cmap=plt.cm.nipy_spectral, edgecolor="k")
plt.show()
# plt.savefig('pca_plot.png')
sys.exit()

df = pd.DataFrame()
trainData = []
alignmentCutoff = 200
labels = []
for acc in kinases:
    for domainNum in kinases[acc].domains:
        if len(kinases[acc].domains[domainNum]) >= alignmentCutoff:
            labels.append(kinases[acc].group)
            '''
            data, AA = oneHotEncoding(acc, domainNum)
            df = pd.DataFrame(data, columns=AA)
            df['HMM'] = fetchPkinase(acc, domainNum)
            df = pd.concat([df, fetchSeqFeat(acc, domainNum), fetchStrucFeat(acc, domainNum)], axis=1)
            #print (acc, kinases[acc].hmm)
            df['Group'] = kinases[acc].group
            trainData.append(df.to_numpy())
            '''
        #break

print (df)
trainData = np.array(trainData)
sys.exit()
'''
np.stack(trainData, axis=0)
print (trainData.shape)
np.save('trainData.npy', trainData)
numSamples, numRows, numCols = trainData.shape
'''
#for data in trainData:
#    print (data.shape)


# In[31]:


trainData = np.load('trainData.npy')
trainData = trainData[:, :, :28]
numSamples, numRows, numCols = trainData.shape
print (trainData.shape)


# In[55]:


## Create Sampling layer
class Sampling(layers.Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

## Built encoder
latent_dim = 50

encoder_inputs = keras.Input(shape=(numRows, numCols, 1))
x = layers.Conv2D(20, 3, activation="relu", strides=2, padding="same")(encoder_inputs)
x = layers.BatchNormalization()(x)
x = layers.Conv2D(50, 3, activation="relu", strides=2, padding="same")(x)
x = layers.BatchNormalization()(x)
x = layers.Flatten()(x)
x = layers.Dense(20, activation="relu")(x)
x = layers.BatchNormalization()(x)
z_mean = layers.Dense(latent_dim, name="z_mean")(x)
z_log_var = layers.Dense(latent_dim, name="z_log_var")(x)
z = Sampling()([z_mean, z_log_var])
encoder = keras.Model(encoder_inputs, [z_mean, z_log_var, z], name="encoder")
encoder.summary()


## Built decoder
latent_inputs = keras.Input(shape=(latent_dim,))
x = layers.Dense(66 * 7 * 50, activation="relu")(latent_inputs)
x = layers.BatchNormalization()(x)
x = layers.Reshape((66, 7, 50))(x)
x = layers.Conv2DTranspose(50, 3, activation="relu", strides=2, padding="same")(x)
x = layers.BatchNormalization()(x)
x = layers.Conv2DTranspose(20, 3, activation="relu", strides=2, padding="same")(x)
x = layers.BatchNormalization()(x)
#decoder_outputs = layers.Cropping2D(cropping=((0, 0), (2, 2)))(x)
#decoder_outputs = layers.Conv2DTranspose(1, 3, activation="sigmoid", padding="same")(decoder_outputs)
decoder_outputs = layers.Conv2DTranspose(1, 3, activation="sigmoid", padding="same")(x)
decoder = keras.Model(latent_inputs, decoder_outputs, name="decoder")
decoder.summary()


# In[56]:


## Define VAE
class VAE(keras.Model):
    def __init__(self, encoder, decoder, **kwargs):
        super(VAE, self).__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self.total_loss_tracker = keras.metrics.Mean(name="total_loss")
        self.reconstruction_loss_tracker = keras.metrics.Mean(
            name="reconstruction_loss"
        )
        self.kl_loss_tracker = keras.metrics.Mean(name="kl_loss")

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.kl_loss_tracker,
        ]

    def train_step(self, data):
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = self.encoder(data)
            reconstruction = self.decoder(z)
            '''
            reconstruction_loss = tf.reduce_mean(
                tf.reduce_sum(
                    keras.losses.binary_crossentropy(data, reconstruction), axis=(1, 2)
                )
            )
            '''
            reconstruction_loss = tf.reduce_mean(
                tf.reduce_sum(
                    tf.square(
                        tf.subtract(data, reconstruction)
                    )
                )
            )
            kl_loss = -0.5 * (1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
            kl_loss = tf.reduce_mean(tf.reduce_sum(kl_loss, axis=1))
            total_loss = reconstruction_loss + kl_loss
        grads = tape.gradient(total_loss, self.trainable_weights)
        self.optimizer.apply_gradients(zip(grads, self.trainable_weights))
        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.kl_loss_tracker.update_state(kl_loss)
        return {
            "loss": self.total_loss_tracker.result(),
            "reconstruction_loss": self.reconstruction_loss_tracker.result(),
            "kl_loss": self.kl_loss_tracker.result(),
        }

print (trainData.shape)
trainData = trainData.reshape(numSamples, numRows, numCols, 1)
vae = VAE(encoder, decoder)
vae.compile(optimizer=keras.optimizers.Adam())
vae.fit(trainData, epochs=50, batch_size=30)


# In[57]:


## dic = {}
colors = ['red', 'yellow', 'green', 'blue', 'grey', 'orange', 'blue', 'cyan', 'brown', 'gold', 'pink']
setLabels = list(set(labels))
for label, color in zip(setLabels, colors):
    dic[label] = color

labelColors = []
for label in labels:
    labelColors.append(dic[label])


def plot_label_clusters(vae, data, labels=None):
    # display a 2D plot of the digit classes in the latent space
    z_mean, _, _ = vae.encoder.predict(data)
    plt.figure(figsize=(12, 10))
    plt.scatter(z_mean[:, 0], z_mean[:, 1], c=labelColors)
    plt.colorbar()
    plt.xlabel("z[0]")
    plt.ylabel("z[1]")
    plt.show()


#(x_train, y_train), _ = keras.datasets.mnist.load_data()
#x_train = np.expand_dims(x_train, -1).astype("float32") / 255

plot_label_clusters(vae, trainData, labelColors)


# In[ ]:
