# -------- THE FOLLOWING PART IS TO GENERATE THE RULES ---------
# Scroll down to line 47 ish unless you want to redo this part
# Run all of these libraries
library(rmcfs)
library(devtools)
library(R.ROSETTA)
library(tidyverse)


# Reading data (must be in same folder)
data <- read.table("Project4.csv", header=TRUE)

# Drops columns 2 through 7 and keeps everything else + remove column 1258 in the original (since 6 columns were removed before it).
genes1 <- data[-c(2:7)]
genes <- genes1[c(-1252)]
table(genes$ethnicity) # display how many rows that are connected to the decisions

# Monte Carlo Feature Selection and Interdependency Discovery. 
# For ranking features by how useful they are for predicting ethnicity.
# This takes ca 10 min
result <- mcfs(ethnicity ~ ., genes, projections = 1500, projectionSize = 0.1, splits = 5, splitSetSize = 500, cutoffPermutations = 6, threadsNumber = 8)

# Render results and save them
plot(result, type = "distances") # mcfs convergence graph
saveRDS(result, file = "mcfs_result_project4.rds")
result <- readRDS("mcfs_result_project4.rds")

# Create vector of statistically significant genes
sig <- result$RI$attribute[1:result$cutoff_value]

# Modify dataframe to move decision last to just keep the significant ones
sig_data <- genes[, c("object_id", sig, "ethnicity")]

# Run rosetta (exclude the first column (object id))
ros <- rosetta(sig_data[, -1], roc = TRUE) #accuracy is better here without that enhanced something that they talked about in the lecture
saveRDS(ros, file = "rosetta_results.rds")
write.csv2(ros$main, "rules.csv", row.names = FALSE)
rules <- ros$main
saveRDS(rules, file = "rules_project4.rds")
fwrite(rules,'rules_input_to_cytoscape_function.csv', sep = ";")

# ------- DONE WITH GENERATING RULES --------
# All necessary files are in this folder anyway.



# ------ THE FOLLOWING PART USES THE RULES FROM ROSETTA AND CONNECTS TO CYTOSCAPE TO CREATE 3D VIEW --------
library(RCy3)
library(devtools)


# The code for connecting to github (and that supports the 3D visualisation)
devtools::install_github("tildawiklund/Cytoscape-3D-visualisation")
library(csVisuNet)

# Open cytoscape app

# Visunet visualization
cytoscapePing() # to see if Cytoscape is open

# Switch to Cy3D 3D view
commandsRun("cy3d set renderer rendererID=org.baderlab.cy3d") 

# Generate info for Cytoscape + connect
vis <- visunetcyto(rules[1:100, ]) 


# This part is optional to check on the node info used in Cytoscape
library(data.table)
fwrite(vis$all$nodes,'see_node_info_for_cytoscape.csv', sep = ";")



# To see what version of Cytoscape you are on
# cytoscapeVersionInfo ()

