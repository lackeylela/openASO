import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import tensorflow as tf
import keras

def parse_args():
    parser = argparse.ArgumentParser(description = "Run regression on tab file")
    required = parser.add_argument_group("Required arguments to run")
    required.add_argument("input", type=argparse.FileType('r'))
    args = parser.parse_args()
    return args


def to_one_hot(sequence, base2idx, max_sequence_len):
    """ Create one hot encoding numpy array of shape=(size of vocab, length of seq)
    from an iterable sequence and vocab to index mapping
    """
    one_hot_position = np.zeros(shape=(4, max_sequence_len), dtype=np.int32)
    for position, base in enumerate(sequence):
        one_hot_position[base2idx[base], position] = 1

    return one_hot_position


def get_input_file():
    """ Create X numpy array of shape=(# of samples, 4, length of sequence) and
    Y numpy array of shape=(# of samples,) from the text input.
    """
    base2idx = {"A":0, "T":1, "G":2, "C":3}

    with open(args.input.name) as datafile:

        lines = [line for line in datafile][1:] #remove heading line

        entry_count = len(lines)
        max_sequence_len = max([len(line) for line in lines])

        #TODO
        # Change one hot encoding into depth 4 instead
        X = np.zeros(shape=(entry_count, 4, max_sequence_len)) # 4 bases
        Y = np.zeros(shape=(entry_count, ))

        for index, line in enumerate(lines):
            columns = line.split()

            sequence = columns[-2]
            effectiveness = columns[-1]

            one_hot_position = to_one_hot(sequence, base2idx, max_sequence_len)

            X[index, :] = one_hot_position
            Y[index] = effectiveness

        return X, Y


# https://keras.io/examples/vision/grad_cam/
def make_gradcam_heatmap(
    img_array, model, last_conv_layer_name, classifier_layer_names
):
    # First, we create a model that maps the input image to the activations
    # of the last conv layer
    last_conv_layer = model.get_layer(last_conv_layer_name)
    last_conv_layer_model = keras.Model(model.inputs, last_conv_layer.output)

    # Second, we create a model that maps the activations of the last conv
    # layer to the final class predictions
    classifier_input = keras.Input(shape=(last_conv_layer.output.shape.as_list()[1:]))
    x = classifier_input
    for layer_name in classifier_layer_names:
        x = model.get_layer(layer_name)(x)
    classifier_model = keras.Model(classifier_input, x)

    # Then, we compute the gradient of the top predicted class for our input image
    # with respect to the activations of the last conv layer
    with tf.GradientTape() as tape:
        # Compute activations of the last conv layer and make the tape watch it
        last_conv_layer_output = last_conv_layer_model(img_array)
        tape.watch(last_conv_layer_output)
        # Compute class predictions
        #TODO:
        # CHANGE THIS INTO REGRESSION
        preds = classifier_model(last_conv_layer_output)
        top_pred_index = tf.argmax(preds[0])
        top_class_channel = preds[:, top_pred_index]

    # This is the gradient of the top predicted class with regard to
    # the output feature map of the last conv layer
    grads = tape.gradient(top_class_channel, last_conv_layer_output)

    # This is a vector where each entry is the mean intensity of the gradient
    # over a specific feature map channel
    pooled_grads = tf.reduce_mean(grads, axis=(0, 1, 2))

    # We multiply each channel in the feature map array
    # by "how important this channel is" with regard to the top predicted class
    last_conv_layer_output = tf.keras.backend.eval(last_conv_layer_output)[0]
    pooled_grads = tf.keras.backend.eval(pooled_grads)
    for i in range(pooled_grads.shape[-1]):
        last_conv_layer_output[:, :, i] *= pooled_grads[i]

    # The channel-wise mean of the resulting feature map
    # is our heatmap of class activation
    heatmap = np.mean(last_conv_layer_output, axis=-1)

    # For visualization purpose, we will also normalize the heatmap between 0 & 1
    heatmap = np.maximum(heatmap, 0) / np.max(heatmap)
    return heatmap


def main():

    X, Y = get_input_file()
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

    # Add a depth (channel) dimension
    X_train =  np.expand_dims(X_train, axis=-1)
    X_test =  np.expand_dims(X_test, axis=-1)

    #This is what a example image conversion looks like
    # plt.figure()
    # plt.imshow(X[0])
    # plt.show()

    shape = X_train[0].shape

    inputs = keras.Input(shape=shape)
    x = keras.layers.Conv2D(3, kernel_size=(4, 2), padding='same', name="conv1")(inputs)
    x = keras.layers.Conv2D(5, kernel_size=(1, 1), padding="same", name="conv2")(x)
    x = keras.layers.MaxPooling2D(pool_size=(2, 2), name="max_pool")(x)
    x = keras.layers.Flatten(name="flatten")(x)
    x = keras.layers.Dense(32, activation='relu', name="dense1")(x)
    outputs = keras.layers.Dense(1, activation="sigmoid", name="dense2")(x)

    model = keras.Model(inputs=inputs, outputs=outputs)
    model.compile(
        loss="binary_crossentropy",
        optimizer='adam',
        metrics=[tf.keras.metrics.MeanSquaredError()]
    )
    model.fit(
        X_train,
        Y_train,
        batch_size=32,
        epochs=100,
        validation_data=(X_test, Y_test)
    )

    last_conv_layer_name = "conv2"
    classifier_layer_names = [
        "max_pool",
        # "dropout1",
        "flatten",
        "dense1",
        # "dropout2",
        "dense2"
    ]

    examples = tf.convert_to_tensor(X_test, dtype=tf.float32)
    heatmap = make_gradcam_heatmap(
        examples,
        model,
        last_conv_layer_name,
        classifier_layer_names
    )

    plt.matshow(heatmap[:, :22])
    plt.show()


if __name__ == "__main__":
    args = parse_args()
    if args.input:
        main()
    else:
        print("Enter input file")
