import numpy as np
import gzip
from tqdm import tqdm
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Dense, Flatten, Reshape, Input,Conv2D,MaxPooling2D,UpSampling2D
import pandas as pd
import plotly.express as px
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans

class Kinase:
    '''A class to store information about a kinase'''
    def __init__(self, acc):
        self.acc = acc
        self.sequence = None
        self.fasta = ''
        self.aln2seq = {}
        self.seq2aln = {}
        self.seq2pfam = {}
        self.pfam2seq = {}
        self.mutations = {'A':[], 'D':[], 'R':[]}
    
    def find_fasta(self):
        # print (self.acc)
        for line in gzip.open('../../KA/UniProtFasta2/'+self.acc+'.fasta.gz', 'rt'):
            if line[0] == '>':
                continue
            self.fasta += line.strip()

def get_train_set2():
    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    df = pd.read_csv('../trainDataFromTrimmedAln.tsv.gz', sep = '\t')

    # Enable this if you want to plot only train data
    # df = df[df.Dataset == 'train']
    # print (df.to_numpy().tolist())
    # sys.exit()

    df['Dataset'] = df['Dataset'].replace(to_replace='train', value=0.025, regex=True)
    df['Dataset'] = df['Dataset'].replace(to_replace='test', value=0.3, regex=True)
    # exclude columns
    # df = df.loc[:, ~df.columns.isin(['allHomologs','exclParalogs','specParalogs','orthologs', 'bpso','bpsh'])]
    df = df.loc[:, ~df.columns.isin([
                                # 'allHomologs',
                                'exclParalogs',
                                'specParalogs',
                                'orthologs'
                                'bpso',
                                'bpsh'
                                ])]
    # exclude columns to make the data matrix
    original_df = df.copy()
    columns_to_exclude = [
                        'Acc',
                        'Mutation',
                        'Gene',
                        'Dataset',
                        'hmmPos',
                        'hmmSS',
                        #   'A_known',
                        #   'D_known',
                        #   'R_known',
                        #   'Phosphomimic',
                        #   'hmmScoreWT',
                        #   'hmmScoreMUT',
                        #   'hmmScoreDiff'
                        ]
    # for aa in AA:
    #     columns_to_exclude.append(aa+'_WT')
    #     columns_to_exclude.append(aa+'_MUT')
    pfam_ptm_cols = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
    for i in range(-5,6):
        if i in [-1, 0, 1]: continue
        for col in pfam_ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

    # ptm_cols = ['p', 'ac', 'me', 'gl', 'm1', 'm2', 'm3', 'sm', 'ub']
    # for i in range(-5,6):
    #     if i in [-1, 0, 1]: continue
    #     for col in pfam_ptm_cols:
    #         columns_to_exclude.append(col.split('_')[0]+'_'+str(i))

    adr_cols = ['A', 'D', 'R']
    for i in range(-5, 6):
        if i in [-1, 0, 1]: continue
        for col in adr_cols:
            columns_to_exclude.append(col+'_'+str(i))

    df = df.loc[:, ~df.columns.isin(columns_to_exclude)]

    discrete_map={
                "A": 1,
                "D": 2,
                "R": 3,
                "N": 4,
                "AR": 5,
                "TBD": 0,
                "Activating": 6,
                "Inconclusive": 7
                }
    color_discrete_map={
                        "A": "green",
                        "D": "red",
                        "R": "blue",
                        "N": "cyan",
                        "AR": "violet",
                        "TBD": "yellow",
                        "Activating": "lightgreen",
                        "Inconclusive": "grey"
                        }
    y_train = []
    y_train_color = []
    y_train_radius = []
    for y in df['MUT_TYPE'].to_numpy():
        y_train.append(discrete_map[y])
        y_train_color.append(color_discrete_map[y])
        if y in ['TBD', 'Activating', 'Inconclusive']:
            y_train_radius.append(100)
        else:
            y_train_radius.append(10)
    y_train = np.array(y_train)
    y_train_color = np.array(y_train_color)
    print (y_train)

    data = df.to_numpy()[:,:-1]
    scaler = MinMaxScaler()
    scaler.fit(data)
    data = scaler.transform(data)
    print (df)
    x_train = []
    for row in data:
        row = [float(item) for item in row[:-1]]
        data = [row]
        x_train.append(data)
    x_train = np.array(x_train)
    
    return x_train, y_train, y_train_color, y_train_radius

        
x_train, y_train, y_train_color, y_train_radius = get_train_set2()
x_test, y_test, y_test_color, y_test_radius = x_train, y_train, y_train_color, y_train_radius
# exit()

# Load the MNIST data set
# (x_train, y_train), (x_test, y_test) = tf.keras.datasets.mnist.load_data()
print (x_train[0], y_train[0])
print (x_train.shape)
# exit()

# Normalize pixel values to [0., 1.]
# x_train = x_train / 255.
# x_test = x_test / 255.

latent_dim = 16

# Images are 28 by 28
img_shape = (x_train.shape[1], x_train.shape[2])

encoder = Sequential([
    Flatten(input_shape=img_shape),
    # Dense(256, activation='relu'),
    Dense(128, activation='relu'),
    Dense(64, activation='relu'),
    Dense(32, activation='relu'),
    Dense(latent_dim, name='encoder_output')
])

decoder = Sequential([
    Dense(32, activation='sigmoid', input_shape=(latent_dim,)),
    Dense(64, activation='sigmoid'),
    Dense(128, activation='sigmoid'),
    # Dense(256, activation='sigmoid'),
    Dense(img_shape[0] * img_shape[1], activation='sigmoid'),
    Reshape(img_shape)
])

'''
encoder = Sequential([
    #input = 28 x 28 x 1 (wide and thin)
    Conv2D(32, (3, 3), activation='relu', padding='same'), #28 x 28 x 32
    MaxPooling2D(pool_size=(2, 2)), #14 x 14 x 32
    Conv2D(64, (3, 3), activation='relu', padding='same'), #14 x 14 x 64,
    MaxPooling2D(pool_size=(2, 2)), #7 x 7 x 64,
    Conv2D(128, (3, 3), activation='relu', padding='same') #7 x 7 x 128 (small and thick)
])

decoder = Sequential([
    Conv2D(128, (3, 3), activation='relu', padding='same'), #7 x 7 x 128
    UpSampling2D((2,2)), # 14 x 14 x 128
    Conv2D(64, (3, 3), activation='relu', padding='same'), # 14 x 14 x 64
    UpSampling2D((2,2)), # 28 x 28 x 64
    Conv2D(1, (3, 3), activation='sigmoid', padding='same') #28 x 28 x 1
])
'''

class TestEncoder(tf.keras.callbacks.Callback):
    def __init__(self, x_test, y_test, y_test_color, y_test_radius):
        super(TestEncoder, self).__init__()
        self.x_test = x_test
        self.y_test = y_test
        self.y_test_color = y_test_color
        self.y_test_radius = y_test_radius
        self.current_epoch = 0

    def on_epoch_begin(self, epoch, logs={}):
        if self.current_epoch%25 == 0: print (self.current_epoch)
        self.current_epoch = self.current_epoch + 1
        encoder_model = Model(inputs=self.model.input,
                              outputs=self.model.get_layer('encoder_output').output)
        encoder_output = encoder_model(self.x_test)
        # plt.subplot(25, 4, self.current_epoch)
        if self.current_epoch == 100:
            color_discrete_map={
                        "A": "green",
                        "D": "red",
                        "R": "blue",
                        "N": "cyan",
                        "AR": "violet",
                        "TBD": "yellow",
                        "Activating": "lightgreen",
                        "Inconclusive": "grey"
                        }
            plt.subplot(1, 1, 1)
            plt.scatter(encoder_output[:, 0],
                        encoder_output[:, 1], s=self.y_test_radius[0:x_test.shape[0]], alpha=0.8,
                        cmap=color_discrete_map, c=self.y_test_color[0:x_test.shape[0]])
            # plt.xlim(-1, 1)
            # plt.ylim(-1, 1)
            plt.xlabel('Latent Dimension 1')
            plt.ylabel('Latent Dimension 2')
            X = []
            for ld1, ld2 in zip(encoder_output[:, 0].numpy(), encoder_output[:, 1]):
                row =[]
                row.append(ld1)
                row.append(ld2)
                X.append(np.array(row))
            X = np.array(X)
            kmeans = KMeans(n_clusters=5, random_state=0, n_init='auto').fit(X)
            print (list(kmeans.cluster_centers_))


autoencoder = Model(inputs=encoder.input, outputs=decoder(encoder.output))
autoencoder.compile(loss='binary_crossentropy', optimizer='adam')

plt.figure(figsize=(15,15))
model_history = autoencoder.fit(x_train, x_train, epochs=100, batch_size=32, verbose=0,
                                callbacks=[TestEncoder(x_test, y_test, y_test_color, y_test_radius)])
plt.show()

plt.plot(model_history.history["loss"])
plt.title("Loss vs. Epoch")
plt.ylabel("Loss")
plt.xlabel("Epoch")
plt.grid(True)
plt.show()
exit()

# (x_train, _), (x_test, _) = keras.datasets.mnist.load_data()
# print (x_train)
x_train = get_train_set()
# exit()
# mnist_digits = np.concatenate([x_train, x_test], axis=0)
# mnist_digits = np.expand_dims(mnist_digits, -1).astype("float32") / 255

vae = VAE(encoder, decoder)
vae.compile(optimizer=keras.optimizers.Adam())
# vae.fit(mnist_digits, epochs=30, batch_size=128)
vae.fit(x_train, epochs=30, batch_size=128)

import matplotlib.pyplot as plt


def plot_latent_space(vae, n=30, figsize=15):
    # display a n*n 2D manifold of digits
    digit_size = 28
    scale = 1.0
    figure = np.zeros((digit_size * n, digit_size * n))
    # linearly spaced coordinates corresponding to the 2D plot
    # of digit classes in the latent space
    grid_x = np.linspace(-scale, scale, n)
    grid_y = np.linspace(-scale, scale, n)[::-1]

    for i, yi in enumerate(grid_y):
        for j, xi in enumerate(grid_x):
            z_sample = np.array([[xi, yi]])
            x_decoded = vae.decoder.predict(z_sample)
            digit = x_decoded[0].reshape(digit_size, digit_size)
            figure[
                i * digit_size : (i + 1) * digit_size,
                j * digit_size : (j + 1) * digit_size,
            ] = digit

    plt.figure(figsize=(figsize, figsize))
    start_range = digit_size // 2
    end_range = n * digit_size + start_range
    pixel_range = np.arange(start_range, end_range, digit_size)
    sample_range_x = np.round(grid_x, 1)
    sample_range_y = np.round(grid_y, 1)
    plt.xticks(pixel_range, sample_range_x)
    plt.yticks(pixel_range, sample_range_y)
    plt.xlabel("z[0]")
    plt.ylabel("z[1]")
    plt.imshow(figure, cmap="Greys_r")
    plt.show()

plot_latent_space(vae)