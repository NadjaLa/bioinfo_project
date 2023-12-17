#Naive Bayes aim + visualizations

install.packages("e1071")
library(e1071)

# combine data frames and make a factor variable for the class column
combined_df <- rbind(Bird_worksheet, mammals_worksheet)
combined_df$Class <- as.factor(combined_df$Class)

# get features from the sequence column
combined_df$QNE <- grepl("QNE", combined_df$sequence)
combined_df$KNE <- grepl("KNE", combined_df$sequence)
combined_df$LGLF <- grepl("LGLF", combined_df$sequence)
combined_df$LAIF <- grepl("LAIF", combined_df$sequence)
combined_df$LVANFR <- grepl("LVANFR", combined_df$sequence)
combined_df$MVEDFR <- grepl("MVEDFR", combined_df$sequence)
combined_df$SKIK <- grepl("SKIK", combined_df$sequence)
combined_df$ANF <- grepl("ANF", combined_df$sequence)
combined_df$ENF <- grepl("ENF", combined_df$sequence)
combined_df$QDKD <- grepl("QDKD", combined_df$sequence)
combined_df$GYAQ <- grepl("GYAQ", combined_df$sequence)
combined_df$GYPE <- grepl("GYPE", combined_df$sequence)
combined_df$ELHD <- grepl("ELHD", combined_df$sequence)
combined_df$WAIL <- grepl("WAIL", combined_df$sequence)
combined_df$WSVL <- grepl("WSVL", combined_df$sequence)

# split into training/testing
set.seed(123) 
sample_index <- sample(1:nrow(combined_df), 0.2 * nrow(combined_df))
train_data <- combined_df[sample_index, ]
test_data <- combined_df[-sample_index, ]

# build Naive Bayes model
features <- c("QNE", "KNE", "LGLF", "LAIF", "LVANFR", "MVEDFR", "SKIK", "ANF", "ENF", "QDKD", "GYAQ", "GYPE", "ELHD", "WAIL", "WSVL")
nb_model <- naiveBayes(Class ~ ., data = train_data[, c(features, "Class")])
predictions <- predict(nb_model, test_data)

# evaluate accuracy
confusion_matrix <- table(predictions, test_data$Class)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Accuracy:", accuracy))

# make predictions for full data set
combo_predictions <- predict(nb_model, newdata = combined_df)

# evaluate accuracy
combo_confusion_matrix <- table(combo_predictions, combined_df$Class)
combo_accuracy <- sum(diag(combo_confusion_matrix)) / sum(combo_confusion_matrix)
print(paste("Accuracy on combo:", combo_accuracy))

### ATTEMPS TO VISUALIZE

# heat map
library(ggplot2)
library(reshape)

confusion_df <- as.data.frame(as.table(combo_confusion_matrix))

ggplot(confusion_df, aes(x = combo_predictions, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Confusion Matrix Heatmap", x = "Prediction", y = "Reference")

# ROC
library(pROC)
install.packages("PRROC")
library(PRROC)

roc_curve <- roc(test_data$Class, as.numeric(predictions == "positive_class"))
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)

# PR curves for model and full data set
pr_values <- pr.curve(test_data$Class, as.numeric(predictions == "positive_class"), curve = TRUE)
plot(pr_values, main = "Precision-Recall Curve", col = "green", lwd = 2)

pr_values2 <- pr.curve(combo_predictions, as.numeric(predictions == "positive_class"), curve = TRUE)
plot(pr_values2, main = "Precision-Recall Curve 2", col = "green", lwd = 2)





