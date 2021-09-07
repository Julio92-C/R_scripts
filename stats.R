
# Run stats analyis of amra data frame

# Basic stats analysis ###########
basic_stats <- function(data_sets) {                         
  cat("Basic stats analysis!")
	sink('./test/Stats_Summary.txt')

	cat("=======================================================================================================\n")
	cat("Stats Summary of the AMR data sets\n")
	cat("=======================================================================================================\n")
	m <- summary(data_sets)
  print(m)
  cat("=======================================================================================================\n")
  
  ## Pearson correlation between two variables
  cat("Pearson correlation between Accuracy and Identity\n")
  pear_corr <- cor(Accuracy, Identity, method = "pearson")
  print(pear_corr)
  cat("=======================================================================================================\n")
	
  # Linear regression model Accuracy vs Identity
  cat("Linear regression model Accuracy vs Identity\n")
  model <- lm(Accuracy~Identity, data = data_sets)
  model_sum <- summary(model)
  print(model_sum)
  cat("=======================================================================================================\n")
  
  # Comparison - T-test
  cat("T-test are used when the two set of population data are normally distributed\n")
  cat("data1(Accuracy, Model_Type == rRNA mutation model)\n")
  cat("data2(Accuracy, Model_Type == protein homolog model)\n")
  data1 <- subset(Accuracy, Model_Type == "rRNA mutation model")
  data2 <- subset(Accuracy, Model_Type == "protein homolog model")
  T_test <- t.test(data1, data2)
  print(T_test)
  cat("=======================================================================================================\n")
  
  # Stop writing to the file
	sink()	
}


