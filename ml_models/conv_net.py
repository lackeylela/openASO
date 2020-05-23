from keras.preprocessing import sequence
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.layers import InputLayer
from keras.layers import Conv1D, GlobalMaxPooling1D
from tensorflow import keras
import tensorflow as tf
import numpy as np

class OneDConvNet:
    def __init__(self):
        strides = 1
        hidden_dims = 30

        model = Sequential()
        model.add(Conv1D(filters=3,
                         kernel_size=3,
                         padding='valid',
                         activation='relu',
                         strides=1,
                         input_shape=(88, 1)))

        model.add(Conv1D(filters=5,
                         kernel_size=3,
                         padding='valid',
                         activation='relu',
                         strides=2))

        model.add(Conv1D(filters=10,
                         kernel_size=3,
                         padding='valid',
                         activation='relu',
                         strides=2))

        model.add(GlobalMaxPooling1D())

        model.add(Dense(hidden_dims))
        # model.add(Dropout(0.2))
        model.add(Activation('relu'))

        model.add(Dense(1))
        model.add(Activation('sigmoid'))

        self.model = model


    def fit(self, X_train, X_test, Y_train, Y_test, epochs=1000, batch_size=32):

        X_train =  np.expand_dims(X_train, axis=-1)
        X_test =  np.expand_dims(X_test, axis=-1)
        Y_train =  np.expand_dims(Y_train, axis=-1)
        Y_test =  np.expand_dims(Y_test, axis=-1)

        self.model.compile(loss="binary_crossentropy",
                      optimizer='adam',
                      metrics=[tf.keras.metrics.MeanSquaredError()])

        self.model.fit(X_train, Y_train,
                  batch_size=batch_size,
                  epochs=epochs,
                  validation_data=(X_test, Y_test))

        return self.model
