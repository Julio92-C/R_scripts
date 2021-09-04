
#!/usr/bin/env Rscript --vanilla

# Script to access the stats of AMR profile

# Author: Julio Cesar Ortega Cambara 

# Date: 30th August, 2021

# Version: v1.0.0

# Load the packages #####
library(ggplot2)
library(ggpubr)
library(argparser)
library(plotly)
library(manipulateWidget)

# Use of argparser in the script
# 1. Firstly we set up a new argument parser, and give it a name:
p <- arg_parser("AMR SUMMARY STATS")

# 2. Then we can add the input_file argument. This argument is mandatory, so we haven’t specified a default

p <- add_argument(p, "input_file", help="Path to input_file, it shoul be a comma separated csv file from  Metrichor’s Antimicrobial Resistance Mapping Application (ARMA)")


# 3. Flag argument. This is a optional argument, that can be used to turn functionality within our script on and off.
p <- add_argument(p, "--filter", help="Filter by Accuracy, Identity, Coverage, Sequence Length", flag=TRUE)
p <- add_argument(p, "--output_dir", help="Put the results to a specific output directory", flag=TRUE)

# 4. parse the arguments:
argv <- parse_args(p)


# Load the data sets ###############################
data_sets <- read.csv(argv$input_file)

# Clean up Data sets ##############################################
data_sets_cleaned <- data_sets[c(-1,-2,-3,-4,-5,-10,-11, -13, -15, -18, -19)]

# Reordering the columns in a data frame #####
data_sets_cleaned <- data_sets_cleaned[c(1, 2, 5, 7, 3, 4, 6, 8)]


# Filter data-set by (identity > 80%, Accurracy > 80% and Coverage > 80% ) ############################## 
data_sets_filtered <- data_sets_cleaned[data_sets_cleaned$identity > 80 & data_sets_cleaned$accuracy > 80
                                        & data_sets_cleaned$coverage > 80,]


# Set up  Variables ######################
Accuracy <- data_sets_filtered$accuracy
Identity <- data_sets_filtered$identity
Coverage <- data_sets_filtered$coverage
Seq_length <- data_sets_filtered$seq_len
Model_Type <- factor(data_sets_filtered$model)
Category <-  factor(data_sets_filtered$categories)
AMR <- data_sets_filtered$name
Specie <- data_sets_filtered$taxon


# Basic stats analysis ###########
sink('./test/Stats_Summary.txt')

cat("=======================================================================================================\n")
cat("Stats Summary of the AMR data sets\n")
cat("=======================================================================================================\n")

# Stats summary of AMR data sets
summary(data_sets_filtered)

## Pearson correlation between two variables
cat("=======================================================================================================\n")
cat("Pearson correlation between Accuracy and Identity\n")
cor(Accuracy, Identity, method = "pearson")

# Linear regression model Accuracy vs Identity
cat("=======================================================================================================\n")
cat("Linear regression model Accuracy vs Identity\n")
model <- lm(Accuracy~Identity, data = data_sets_filtered)
summary(model)

# Comparison - T-test
cat("=======================================================================================================\n")
cat("T-test are used when the two set of population data are normally distributed\n")
cat("data1(Accuracy, Model_Type == rRNA mutation model)\n")
cat("data2(Accuracy, Model_Type == protein homolog model)\n")

data1 <- subset(Accuracy, Model_Type == "rRNA mutation model")
data2 <- subset(Accuracy, Model_Type == "protein homolog model")
t.test(data1, data2)


# Stop writing to the file
sink()

# ARM Model Types Frequency
q <- ggplot(data_sets_filtered, aes(x=factor(Model_Type), fill=Model_Type))
q <- q + labs(x="Models", y="Frequency")
q <- q + geom_bar(stat="count", width=0.7, color="black")
q <- q + scale_x_discrete(limits=c("rRNA mutation model", "protein homolog model", 
                                   "protein variant model", "protein wild type model"))
q <- q +  theme_classic()
q <- q + ggtitle("AMR: Model Types Frequency")
q <- q + theme(plot.title = element_text(hjust = 0.5))
q <- q + theme(legend.position="none")
q

# AMR: Identity vs Accuracy
p <- ggplot(data_sets_filtered, aes(x=Accuracy, y=Identity, size=Seq_length, 
                                    shape=Model_Type))
p <- p + geom_point(alpha = 0.5, aes(colour=Model_Type))
p <- p + theme_classic()
p <- p + theme(legend.position = "none")
p <- p + ggtitle("AMR: Identity vs. Accuracy")
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + labs(x="Accurracy", y="Identity")
p <- p + theme(legend.position="none")
p

# Coverage vs sequence length
c <- ggplot(data_sets_filtered, aes(x=Seq_length, y=Coverage, fill=Model_Type,))
c <- c + geom_bin_2d()
c <- c + theme_classic()
c <- c + ggtitle("AMR: Coverage vs sequence length")
c <- c + theme(plot.title = element_text(hjust = 0.5))
c <- c + labs(x="Sequence Length (bp)", y="Coverage (%)")
c <- c + theme(legend.position="none")
c

#  Boxplot() Model type vs Identity
b <- ggplot(data_sets_filtered, aes(x=reorder(Model_Type, Identity, FUN = median), 
                                    y=Identity, fill=factor(Model_Type)))
b <- b + geom_boxplot()
b <- b + stat_boxplot(geom = 'errorbar', width=0.2)
b <- b + theme_classic()
b <- b + ggtitle("AMR: Identity vs Model Type")
b <- b + theme(plot.title = element_text(hjust = 0.5))
b <- b + theme(legend.position="none")
b

#  Boxplot() Model type vs Accuracy
a <- ggplot(data_sets_filtered, aes(x=reorder(Model_Type, Accuracy, FUN = median), y=Accuracy, fill=Model_Type))
a <- a + geom_boxplot()
a <- a + stat_boxplot(geom = 'errorbar', width=0.2)
a <- a + theme_classic()
a <- a + ggtitle("AMR: Accuracy vs Model Type")
a <- a + theme(plot.title = element_text(hjust = 0.5))
a <- a + theme(legend.position="none")
a

# AMR: Category vs Model Type
j <- ggplot(data_sets_filtered, aes(y=factor(Category, labels = "CAT"), x=Model_Type, col=Model_Type, size=Seq_length))
j <- j + labs(x="Model Type", y="AMR Category")
j <- j + geom_point()
j <- j + theme_minimal()
#j <- j + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
j <- j + ggtitle("AMR: Category vs Model Type") 
j <- j + theme(plot.title = element_text(hjust = 0.5))
j <- j + theme(legend.position="none")
j


# Display all of the graph in the same page
ggarrange(c, p, a, b, q, j + rremove("x.text"), 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)


# Save plots to a pdf file
ggsave("./test/amra_stats_summary.pdf", width = 34, height = 38, units = "cm")



# Plot Model type vs Identity with PLOTLY ####
fig1 <- plotly::ggplotly(b)
fig1


# Plot Model type vs Accuracy with PLOTLY ####
fig2 <- plotly::ggplotly(a)
fig2

# Plot ARM Model Types Frequency with PLOTLY ####
fig3 <- plotly::ggplotly(q)
fig3

# AMR: Identity vs Accuracy with PLOTLY ####
fig4 <- plotly::ggplotly(p)
fig4

# AMR: Sequence length vs Coverage ####
fig5 <- plotly::ggplotly(c)
fig5


# Save Graph to HTML file ####
htmlwidgets::saveWidget(fig1, "./test/model_type_accur.html", selfcontained = F)

htmlwidgets::saveWidget(fig2, "./test/model_type_ident.html", selfcontained = F)

htmlwidgets::saveWidget(fig3, "./test/model_type_freq.html", selfcontained = F)

htmlwidgets::saveWidget(fig4, "./test/accur_ident.html", selfcontained = F)

htmlwidgets::saveWidget(fig5, "./test/seq_length_coverage.html", selfcontained = F)

combine_graph <- combineWidgets(fig3, fig2, fig1, fig5)
combine_graph

# Summary of AMR stats ####
htmlwidgets::saveWidget(combine_graph, "./test/amra_summary.html")



# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages
detach("package:datasets", unload = TRUE)
