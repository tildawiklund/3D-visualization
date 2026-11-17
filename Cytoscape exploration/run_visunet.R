library(rmcfs)
library(devtools)
library(R.ROSETTA)
library(tidyverse)
library(csVisuNet)
library(RCy3)

result <- readRDS("mcfs_result_project4.rds")

ros <- readRDS("rosetta_results.rds")

rules <- readRDS("rules_project4.rds")


#Visunet visualization
cytoscapePing()
vis <- visunetcyto(rules[1:100, ])
