rm(list=ls())
library(keras)

mnist <- dataset_mnist()
mnist$train$x <- mnist$train$x/255
mnist$test$x <- mnist$test$x/255

vectorize_sequences <- function(sequences, dimension = c(28, 28)) {
        # Create an all-zero matrix of shape (len(sequences), dimension)
        results <- matrix(0, nrow = length(sequences), ncol = dimension)
        for (i in 1:length(sequences))
                # Sets specific indices of results[i] to 1s
                results[i, sequences[[i]]] <- 1
        results
}

original_model <- keras_model_sequential() %>% 
        layer_dense(units = 16, activation = "relu", input_shape = c(28, 28)) %>% 
        layer_dense(units = 16, activation = "relu") %>% 
        layer_dense(units = 1, activation = "sigmoid")

original_model %>% compile(
        optimizer = "rmsprop",
        loss = "binary_crossentropy",
        metrics = c("accuracy")
)

smaller_model <- keras_model_sequential() %>% 
        layer_dense(units = 4, activation = "relu", input_shape = c(28, 28)) %>% 
        layer_dense(units = 4, activation = "relu") %>% 
        layer_dense(units = 1, activation = "sigmoid")

smaller_model %>% compile(
        optimizer = "rmsprop",
        loss = "binary_crossentropy",
        metrics = c("accuracy")
)




# Our vectorized training data
x_train <- mnist$train$x
# Our vectorized test data
x_test <- mnist$test$x

# Our vectorized labels
y_train <- as.numeric(mnist$train$y)
y_test <- as.numeric(mnist$test$y)
original_hist <- original_model %>% fit(
        x_train, y_train,
        epochs = 20,
        batch_size = 512,
        validation_data = list(x_test, y_test)
)

smaller_model_hist <- smaller_model %>% fit(
        x_train, y_train,
        epochs = 20,
        batch_size = 512,
        validation_data = list(x_test, y_test)
)

library(ggplot2)
library(tidyr)
plot_training_losses <- function(losses) {
        loss_names <- names(losses)
        losses <- as.data.frame(losses)
        losses$epoch <- seq_len(nrow(losses))
        losses %>% 
                gather(model, loss, loss_names[[1]], loss_names[[2]]) %>% 
                ggplot(aes(x = epoch, y = loss, colour = model)) +
                geom_point()
}


plot_training_losses(losses = list(
        original_model = original_hist$metrics$val_loss,
        smaller_model = smaller_model_hist$metrics$val_loss
))

bigger_model <- keras_model_sequential() %>% 
        layer_dense(units = 512, activation = "relu", input_shape = c(28, 28)) %>% 
        layer_dense(units = 512, activation = "relu") %>% 
        layer_dense(units = 1, activation = "sigmoid")

bigger_model %>% compile(
        optimizer = "rmsprop",
        loss = "binary_crossentropy",
        metrics = c('acc')
)

bigger_model_hist <- bigger_model %>% fit(
        x_train, y_train,
        epochs = 20,
        batch_size = 512,
        validation_data = list(x_test, y_test)
)

plot_training_losses(losses = list(
        original_model = original_hist$metrics$val_loss,
        bigger_model = bigger_model_hist$metrics$val_loss
))

plot_training_losses(losses = list(
        original_model = original_hist$metrics$loss,
        bigger_model = bigger_model_hist$metrics$loss
))

l2_model <- keras_model_sequential() %>% 
        layer_dense(units = 16, kernel_regularizer = regularizer_l2(0.001),
                    activation = "relu", input_shape = c(28, 28)) %>% 
        layer_dense(units = 16, kernel_regularizer = regularizer_l2(0.001),
                    activation = "relu") %>% 
        layer_dense(units = 1, activation = "sigmoid")

l2_model %>% compile(
        optimizer = "rmsprop",
        loss = "binary_crossentropy",
        metrics = c("acc")
)

l2_model_hist <- l2_model %>% fit(
        x_train, y_train,
        epochs = 20,
        batch_size = 512,
        validation_data = list(x_test, y_test)
)

plot_training_losses(losses = list(
        original_model = original_hist$metrics$val_loss,
        l2_model = l2_model_hist$metrics$val_loss
))

# L1 regularization
regularizer_l1(0.001)

# L1 and L2 regularization at the same time
regularizer_l1_l2(l1 = 0.001, l2 = 0.001)


dpt_model <- keras_model_sequential() %>% 
        layer_dense(units = 16, activation = "relu", input_shape = c(28, 28)) %>% 
        layer_dropout(rate = 0.5) %>% 
        layer_dense(units = 16, activation = "relu") %>% 
        layer_dropout(rate = 0.5) %>% 
        layer_dense(units = 1, activation = "sigmoid")

dpt_model %>% compile(
        optimizer = "rmsprop",
        loss = "binary_crossentropy",
        metrics = c("acc")
)

dpt_model_hist <- dpt_model %>% fit(
        x_train, y_train,
        epochs = 20,
        batch_size = 512,
        validation_data = list(x_test, y_test)
)

plot_training_losses(losses = list(
        original_model = original_hist$metrics$val_loss,
        dpt_model = dpt_model_hist$metrics$val_loss
))

