####1 Setup
###1.1 Load packages
library(readxl)
library(dplyr)
library(lavaan)
library(car)
library(qgraph)
library(bootnet)
library(networktools)
library(NetworkComparisonTest)
library(igraph)
library(EGAnet)
library(glmnet)
library(psychonetrics)
library(grDevices)
library(psych)
library(tidyverse)
library(showtext)
library(Hmisc)
library(networktools)
library(ggplot2)
library(tidyr)

###1.2.1 Load data
getwd()
setwd("C:\\R\\JYA")

data1<-read_excel("CLPN2_1.xlsx")
data2<-read_excel("CLPN2_2.xlsx")

# Select the required covariates from the T1 data; for example, assume that 'Gender' and 'Age' are available.
covariates_t1 <- data1 %>% select(age, ses) # Replace Gender and Age with the actual covariate names, such as gender.
data_t1<-data1 %>% select(P2_1:WB2_9)## T1 depression
data_t2<-data2 %>% select(P3_1:WB3_9)## T2 depression
# Combine the covariates with the T1 and T2 variable data.
df_with_cov <- cbind(data_t1, data_t2, covariates_t1)

# Set the column names (leave the names of the T1 and T2 variables unchanged).
nt1 <- paste0("T1", colnames(data_t1))
nt2 <- paste0("T2", colnames(data_t2))
nt_cov <- colnames(covariates_t1) # Keep the original covariate names, or add a prefix, such as paste0("Cov_", colnames(covariates_t1)).
colnames(df_with_cov) <- c(nt1, nt2, nt_cov)

# Variable labels (covariates are generally not displayed as network nodes, so the labels may not need to be changed).
labels <- c("P1","P2","P3","P4",
            "N1","N2","N3","N4","N5",
            "WB1","WB2","WB3","WB4","WB5","WB6","WB7","WB8","WB9")

k2 <- 18 # Number of core variables
num_cov <- ncol(covariates_t1) # Number of covariates

# Define the communities first: divide the 27 nodes into three groups of nine.
# Modify the community names based on the actual meaning of the variables.
communities <- c(
  rep("p", 4),
  rep("n", 5),
  rep("W", 9)
)

# Note: The second group in labels contains X5-X9; these may need to be D5-D9.
# If these labels are duplicated, make them unique first; otherwise, the bridge metric plot will be confusing.
names(communities) <- labels

tesds.fun <- function(df){

  adjmat2 <- matrix(0, k2, k2) # Initialize the adjacency matrix with only the core variables.
  for(i in 1:k2){
    # Predictors include all T1 core variables and all covariates.
    predictors <- as.matrix(df[, c(1:k2, (2*k2+1):(2*k2+num_cov))]) 
    # The outcome variable is the ith T2 variable.
    response_var <- df[, k2 + i] 
    
    lassoreg <- cv.glmnet(predictors, response_var,
                          family = "gaussian", alpha = 1,
                          standardize = TRUE, nfolds = 10)
    lambda <- lassoreg$lambda.min
    # When extracting coefficients, retain only the first k2, which correspond to the paths from the T1 core variables to the current T2 variable.
    # For now, exclude the covariate coefficients from the final variable network plot.
    coefs <- as.vector(coef(lassoreg, s = lambda, exact = FALSE))[-1] # Remove the intercept.
    adjmat2[1:k2, i] <- coefs[1:k2] 
  }
  return(adjmat2)
}

threshold <- 0.05
# Estimate the network using the data frame that includes the covariates.
set.seed(1000)
Network_tesds_with_cov <- estimateNetwork(df_with_cov, 
                                          fun = tesds.fun, 
                                           labels = labels, 
                                          directed = TRUE)



# Plot  ################# Network plot
plot(Network_tesds_with_cov, 
     labels = labels, 
     color = "lightblue", 
     edge.label.cex = 1, 
     # edge.labels = TRUE,
     threshold = threshold,
     legend = TRUE)

# Use qgraph to plot the directed cross-lagged network.
# Use the spring layout; the force-directed algorithm automatically places more strongly connected nodes closer together [6](@ref).
cross_lag_matrix <- Network_tesds_with_cov$graph

# 1. Create an edge-label matrix with values rounded to three decimal places.
custom_edge_labels <- matrix(
  sprintf("%.2f", cross_lag_matrix), # Format to three decimal places.
  nrow = nrow(cross_lag_matrix),
  ncol = ncol(cross_lag_matrix)
)

# 2. Set labels for nonsignificant edges to empty strings (using the specified threshold).
custom_edge_labels[abs(cross_lag_matrix) < threshold] <- ""
# Plot using the qgraph package.
pdf("network T23.pdf", width=10, height=10)
qgraph_graph <- qgraph(cross_lag_matrix, threshold=threshold,
                       # edge.labels = custom_edge_labels,
                       edge.label.cex = 1, labels = labels)
dev.off()

####################################################### Edge weights
Boot_edge <- bootnet(  Network_tesds_with_cov,  directed = TRUE,  nCores = 1,  nBoots = 1000,
                       type = "nonparametric",  statistics = "edge")

pdf("Edge weight accuracy 95CI(T23).pdf", width = 6, height = 8)
plot(  Boot_edge,  statistics = "edge",  labels = FALSE,  order = "sample")
dev.off()

####################################################### Centrality
set.seed(1000)
Boot2 <- bootnet(Network_tesds_with_cov, directed = TRUE,nCores = 1, nBoots = 1000,type = "case",
                 statistics = c( "bridgeExpectedInfluence"), communities = communities)

corStability(Boot2,statistics = c( "bridgeExpectedInfluence"))

pdf("CS 中心性指标(T23).pdf", width = 10, height = 7)
plot(  Boot2, statistics = c("bridgeExpectedInfluence"))

dev.off()

########################################################


pdf("EI out and in(T23).pdf", width=6, height=6)
viewt <- centralityPlot(qgraph_graph,labels = labels,scale= c("z-scores"),
                        include=c("InExpectedInfluence","OutExpectedInfluence"))
dev.off()

# Use qgraph to plot the directed cross-lagged network.
cross_lag_matrix <- Network_tesds_with_cov$graph
# ===== Export edge-weight matrix, same style as Table S4/S5 =====
rownames(cross_lag_matrix) <- labels
colnames(cross_lag_matrix) <- labels

edge_weight_matrix <- round(cross_lag_matrix, 3)

edge_weight_table <- data.frame(
  Nodes = rownames(edge_weight_matrix),
  edge_weight_matrix,
  check.names = FALSE
)

print(edge_weight_table)

write.csv(edge_weight_table,"Edge weight matrix T23.csv",  row.names = FALSE)



qgraph_graph <- qgraph(cross_lag_matrix, 
                       # threshold=threshold,
                       edge.label.cex = 1, edge.labels = TRUE, 
                       labels = labels, 
                       DoNotPlot = TRUE)

# Calculate centrality metrics.
centrality_scores <- centrality(qgraph_graph)


############################################################################### Bridge nodes


# Whether to use the thresholded network to calculate bridge metrics
bridge_mat <- cross_lag_matrix
bridge_mat[abs(bridge_mat) < threshold] <- 0
diag(bridge_mat) <- 0
rownames(bridge_mat) <- colnames(bridge_mat) <- labels

# Out bridge expected influence:
# Sum of the signed edge weights from the current node to nodes in other communities
bridge_out <- bridge(
  bridge_mat,
  communities = communities,
  directed = TRUE
)

bridge_EI_1step <- data.frame(
  Node = labels,
  Community = communities,
  Bridge_EI_1step_raw = bridge_out$`Bridge Expected Influence (1-step)`
)

bridge_EI_1step <- bridge_EI_1step %>%
  mutate(
    z_Bridge_EI_1step = as.numeric(scale(Bridge_EI_1step_raw)),
    rank_z_Bridge_EI_1step = rank(-z_Bridge_EI_1step, ties.method = "min")
  )

write.csv(
  bridge_EI_1step,
  "Bridge Expected Influence 1-step z scores(T23).csv",
  row.names = FALSE
)

 ggplot(
  bridge_EI_1step,
  aes(
    x = reorder(Node, z_Bridge_EI_1step),
    y = z_Bridge_EI_1step,
    fill = Community
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  labs(
    x = NULL,
    y = "Bridge Expected Influence (z-score)"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )
