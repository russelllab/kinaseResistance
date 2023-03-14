import numpy as np
import gzip
from tqdm import tqdm
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

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

def get_train_set():
    ## Load the pfam bit scores
    pfam_dic = {}
    for line in open('../../pfam/Pkinase.hmm', 'r'):
        if line.strip() in ['#', '//']:
            continue
        # print (line)
        if line.split()[-2] == '-' and line.split()[-3] == '-':
            pfam_pos = line.split()[0]
            pfam_dic[pfam_pos] = {}
            pfam_bit_scores = [float(num) for num in line.split()[1:21]]
            AA = 'ACDEFGHIKLMNPQRSTVWY'
            for num, pfam_bit_score in enumerate(pfam_bit_scores):
                pfam_dic[pfam_pos][AA[num]] = pfam_bit_score
    
    kinases = {}

    for line in gzip.open('../../data/allKinasesPkinaseMappings.tsv.gz', 'rt'):
        if line[0] == '#':
            continue
        uniprot_position = line.split('\t')[2]
        pfam_position = line.split('\t')[4].strip()
        acc = line.split('\t')[0]
        if '-' in acc:
            continue
        if acc not in kinases:
            kinases[acc] = Kinase(acc)
        kinases[acc].seq2pfam[uniprot_position] = pfam_position
        kinases[acc].pfam2seq[pfam_position] = uniprot_position
        if kinases[acc].fasta == '': kinases[acc].find_fasta()

    '''
    # Resistance Mutations
    for line in gzip.open('../../KA/resistant_mutations_Mar_2023.tsv.gz', 'rt'):
        if line[0] == '#':
            continue
        acc = line.split('\t')[2]
        cosmic_mutation = line.split('\t')[1]
        wtAA = cosmic_mutation[0]
        mutAA = cosmic_mutation[-1]
        uniprot_position = str(line.split('\t')[5].strip())
        uniprot_mutation = wtAA + uniprot_position + mutAA
        mut_type = 'R'
        if uniprot_position not in kinases[acc].mutations[mut_type]:
            kinases[acc].mutations[mut_type].append(uniprot_position)

    # ACT/DEACT mutations
    for line in open('../../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
        if line.split('\t')[0] == 'uniprot_name':
            continue
        acc = line.split('\t')[1]
        mut_type = line.split('\t')[5]
        wtAA = line.split('\t')[2]
        uniprot_position = line.split('\t')[3]
        if uniprot_position not in kinases[acc].mutations[mut_type]:
            kinases[acc].mutations[mut_type].append(uniprot_position)
    '''
    trainData = []
    data = []
    for line in gzip.open('../trainData.tsv.gz', 'rt'):
        if line.split('\t')[0] == 'Acc':
            continue
        row = [float(value) for value in line.split('\t')[6: -1]]
        data.append(np.array(row))
    
    data = np.array(data)
    # trainData.append(data)
    # trainData = np.array(trainData)
    # print (trainData.shape)
    return data
        


class Sampling(layers.Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

latent_dim = 2

'''Encoder'''
encoder_inputs = keras.Input(shape=(503, 154, None))
x = layers.Conv2D(32, 3, activation="relu", strides=2, padding="same")(encoder_inputs)
x = layers.Conv2D(64, 3, activation="relu", strides=2, padding="same")(x)
x = layers.Flatten()(x)
x = layers.Dense(16, activation="relu")(x)
z_mean = layers.Dense(latent_dim, name="z_mean")(x)
z_log_var = layers.Dense(latent_dim, name="z_log_var")(x)
z = Sampling()([z_mean, z_log_var])
encoder = keras.Model(encoder_inputs, [z_mean, z_log_var, z], name="encoder")
encoder.summary()


'''Decoder'''
latent_inputs = keras.Input(shape=(latent_dim,))
x = layers.Dense(1 * 126 * 64, activation="relu")(latent_inputs)
x = layers.Reshape((64, 126, 1))(x)
x = layers.Conv2DTranspose(3, 64, activation="relu", strides=2, padding="same")(x)
x = layers.Conv2DTranspose(2, 32, activation="relu", strides=2, padding="same")(x)
decoder_outputs = layers.Conv2DTranspose(1, 3, activation="sigmoid", padding="same")(x)
decoder = keras.Model(latent_inputs, decoder_outputs, name="decoder")
decoder.summary()

class VAE(keras.Model):
    def __init__(self, encoder, decoder, **kwargs):
        super().__init__(**kwargs)
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
            reconstruction_loss = tf.reduce_mean(
                tf.reduce_sum(
                    keras.losses.binary_crossentropy(data, reconstruction), axis=(1, 2)
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