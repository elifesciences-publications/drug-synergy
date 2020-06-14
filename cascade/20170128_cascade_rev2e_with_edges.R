# Performs various experiments to highlight the synergy
# of T and M in MCF-7 at the drug level.
#
# Authors: Thomas Schaffter (tschaff AT us.ibm.com)
# Last updated: 2017-03-13
#
# 2017-03-13: See "COMMENTED AS REQUESTED BY JENNY" in the middle of the script.

rm(list=ls())

source("R/network/getABUndirectedNetwork.R")
source("R/cascade/getTranscriptionalCascade.R")
source("R/degs/getDEGs.R")

# ====================
# READ NETWORK (only for the list of TFs)

network <- read.csv("data/network/MCF7_Woo_et_al-Li_et_al_genelevel_20160126_JELD.txt",
                    fill = TRUE, sep = "\t", stringsAsFactors = FALSE)
length(TFs <- unique(network$Regulator))

# # ====================
# # UNDIRECTED TRANSCRIPTIONAL NETWORK
# 
# # Extract the transcriptional network
# TFs <- unique(network$Regulator)
# tfNetwork <- network[network$Regulator %in% TFs & network$Target %in% TFs,]
# tfNetwork <- getABUndirectedNetwork(tfNetwork)

# ====================
# UNDIRECTED TRANSCRIPTIONAL NETWORK

# tfNetwork <- readRDS("output/network/network_rds/mcf7_TF_network_ABBA_edgelist.rds")
tfNetwork <- readRDS("output/network/network_rds/mcf7_TF_network_FDR_ABBA_edgelist.rds")

# tfNetwork <- getABUndirectedNetwork(tfNetwork) # required for getTranscriptionalCascade() v1
# TFs <- network

# ====================
# Differentially expressed genes

DEGs <- list(
  name = "limma",
  data = readRDS("data/degs/DEGs_limma.rds"),
  threshold = 1e-18
)
# limmaGenes <- rownames(DEGs$data)

# TFACTs <- list(
#   name = "Fisher",
#   data = readRDS("output/tfacts/TFACTs_negative.rds"),
#   threshold = -1e-20
# )
# DEGs <- TFACTs # use active TFs instead of DEGs

TFACTs <- list(
  name = "Fisher",
  data = -1 * data.matrix(read.table("data/jenny_tfacts/MCF7_aracne_all_final_TFs.tsv",
                         sep = "\t", header = T, row.names = 1, stringsAsFactors = F)),
  threshold = -1e-24 # arbitrary, slightly below 0
)
DEGs <- TFACTs # use active TFs instead of DEGs

# VIPERs <- list(
#   name = "Viper",
#   data = readRDS("data/viper/viper_activity_data.rds"),
#   threshold = -1e-20
# )
# DEGs <- VIPERs # user Viper TFs instead of DEGs

timePoints <- c(0, 3, 6, 9, 12, 24)
replicates <- c("A", "B", "C")
rnaseq <- readRDS("data/rna-seq/counts_log.rds")

drugDegs <- list(
  T = getAllDEGs(DEGs, "T", timePoints),
  M = getAllDEGs(DEGs, "M", timePoints),
  TM = getAllDEGs(DEGs, "TM", timePoints)
)
drugDegs[["all"]] <- Reduce(union, c(drugDegs$T, drugDegs$M, drugDegs$TM))

# Generate the list of new DEGs for every drugs and time points.
newDegs <- list()
drugs <- c("T", "M", "TM")
for (drug in drugs) {
  seen <- NULL
  for (tIndex in 1:length(timePoints)) {
    key <- paste0(drug, "_", timePoints[tIndex])
    if (tIndex == 1) {
      newDegs[[key]] <- getDEGs(DEGs, drug, timePoints[tIndex])
      seen <- newDegs[[key]]
    } else {
      newDegs[[key]] <- setdiff(getDEGs(DEGs, drug, timePoints[tIndex]), seen)
      seen <- union(seen, newDegs[[key]])
    }
  }
}

# Returns the DE code (T-M-TM) of a gene.
getOverallDECode <- function(gene) { #, TFs
  return(paste0(c(gene %in% drugDegs$T, 
                  gene %in% drugDegs$M, 
                  gene %in% drugDegs$TM)*1,
                collapse = ""))
}
allDegsDECode <- sapply(drugDegs$all, getOverallDECode)

# Returns the DE code of a gene for the time point specified.
getTemporalDECode <- function(gene, t) {
  return(paste0(c(gene %in% getDEGs(DEGs, "T", t), 
                  gene %in% getDEGs(DEGs, "M", t), 
                  gene %in% getDEGs(DEGs, "TM", t))*1,
                collapse = ""))
}

# # ====================
# # TFACTs
# 
# 
# 
# TFACTs <- list(
#   name = "Fisher",
#   data = readRDS("output/tfacts/TFACTs_negative.rds"),
#   threshold = -1e-20
# )
# nrow(TFACTs$data)
# sum(rowSums(TFACTs$data) != 0)
# 
# # Keep only TFs from the network
# TFACTs$data <- TFACTs$data[rownames(TFACTs$data) %in% TFs,]
# cat("Number of active TFs that are TFs (expecting all): ", nrow(TFACTs$data), "\n", sep = "")

# ====================
# Static network
#
# - Genes that are not TF in the network are filtered out.

network <- read.csv("data/network/MCF7_Woo_et_al-Li_et_al_genelevel_20160126_JELD.txt",
                    fill = TRUE, sep = "\t", stringsAsFactors = FALSE)
network <- network[network$Spearman_FDR < 0.05,]
length(genes <- union(network$Regulator, network$Target))
length(TFs <- unique(network$Regulator))

# igraph network
gNetwork <- edgelist_to_igraph(network, directed = FALSE)
gNetwork <- igraph::simplify(gNetwork)
gNetwork$name <- "mcf7_network_FDR"

# igraph TF network
gTfNetwork <- delete.vertices(gNetwork, setdiff(genes, TFs))
vcount(gTfNetwork)

# ====================
# Dynamic network

# drug <- "TM"
# timePoints <- c(0, 3, 6, 9, 12, 24)
# 
# # Load dynamic network
# filename <- paste0("output/network/network_rds/mcf7_TF_network_FDR_ABBA_edgelist_", drug, "_dynamic.rds")
# tfNetworks <- readRDS(filename)

# ====================
# Define an empty layer of the cascade.

# The top layer contains the new 101, 011, 111 and the 001 from all
# the previous time points combined together.
topLayer <- list(
  "101" = character(0),
  "011" = character(0),
  "111" = character(0),
  past001 = character(0) # past 001 but not unexplained (== Difference 1 with Rev c ==)
)

# The bottom layer contains the new 001. The 001 are organized into
# different groups depending on how they can be explained.
bottomLayer <- list(
  "001" = character(0),
  groups = NULL
)

# A group of genes and their edges to other set of genes.
emptyEdgeList <- as.data.frame(matrix(NA, nrow = 0, ncol = 2))
group <- list(
  degs = character(0),
  undirEdges = emptyEdgeList
)

# Organization of the new 001 into different groups that
# can explain them.
degs001Layer <- list(
  "111" = group,          # have at least one parent in 111 but not in both 101 and 011
  "101-011" = group,      # have at least one parent in 101 and one in 011, but not in 111
  "101-111-011" = group,  # have at least one parent in 101, 111 and 011
  "001" = NULL,           # see below
  empty = character(0)    # genes that can not be explained by any of the above groups
)

# # Option 1: have at least one parent in one of the three above group
# # Option 2: same than Option 1 but with a recursive search to get
# #           the connected component
# past001 = group,        # New 001 that are connected to the past 001 but without considering
# # the past 001 unexplained (== Difference with Rev c ==)
# # The new 001 here overlap with the ones listed in the previous group


# Small adaptation to separate the edges that come from different sources
degs001Layer$`101-011`$undirEdges <- list(
  source_101 = emptyEdgeList,
  source_011 = emptyEdgeList
)
degs001Layer$`101-111-011`$undirEdges <- list(
  source_101 = emptyEdgeList,
  source_111 = emptyEdgeList,
  source_011 = emptyEdgeList
)
degs001Layer$`001` <- list(
  new = list( # 001 than can be explained by at least one new 001 explained
    # Option 1: only the genes that have a direct connection with
    #           the new 001 that are explained by a motif, the remaining
    #           genes are marked as unexplained
    # Option 2: same than Option 1 but with a recursive search to get
    #           the connected component
    degs = character(0),              
    "111" = list(
      undirEdges = emptyEdgeList
    ),
    "101-111-011" = list(
      undirEdges = emptyEdgeList
    ),
    "101-011" = list(
      undirEdges = emptyEdgeList
    ),
    "001" = list(
      emptyEdgeList # Empty if we don't run the recursive search
    )
  ),
  past = group # 001 that can only be exaplined by past 001
)

bottomLayer$groups <- degs001Layer

# Builds an empty result structure.
resEmpty <- list(
  t = NULL,
  degs = character(0),      # ALL active degs now at at t-1 (not the new degs!)
  newDegs = character(0),
  degsCode = character(0),
  topLayer = topLayer,
  bottomLayer = bottomLayer
)

# ====================
# Helper function

# Returns the neighbors of a set of genes in an igraph network.
getNeighbors <- function(genes, order, gNetwork) {
  vNeighbors <- neighborhood(gNetwork, order, V(gNetwork)[genes])[[1]]
  vNeighbors$name[!vNeighbors$name %in% genes]
}

# Helper function to add edges to an existing data frame (list of edges).
# Example:
# source1 -> target
# source2 -> target
# ...
addSourcesToSingleTargetdEdges <- function(target, sources, df) {
  for (source in sources) {
    df <- rbind(t(as.data.frame(c(source, target))), df)
  }
  rownames(df) <- NULL
  df
}

# ====================
# Identify the cascade
#
# Rev 2e: instead of explaining the new 001 at t, we take 
#         *all* the 001 that are active at t.

recursiveSearch <- FALSE
removeUnexplainedFromPast001 <- TRUE
cascade <- list()
for (tIndex in 1:length(timePoints)) {
  # tIndex <- 6
  
  res <- resEmpty
  
  t <- timePoints[tIndex]
  res$t <- t
  
  # Fills the top layer 101, 111 and 011, and past 001
  # Get all the degs active at time t
  degs_t <- getDEGs(DEGs, "TM", t)
  degs_t <- degs_t[degs_t %in% V(gTfNetwork)$name]
  degs_t_codes <- sapply(degs_t, getTemporalDECode, t)
  res$topLayer$`101` <- names(degs_t_codes[degs_t_codes == "101"])
  res$topLayer$`111` <- names(degs_t_codes[degs_t_codes == "111"])
  res$topLayer$`011` <- names(degs_t_codes[degs_t_codes == "011"])
  
  # # Get all the degs active at time t-1
  # t_prev <- NULL
  # if (length(cascade) > 1) {
  #   layer_prev <- cascade[[length(cascade)]]
  #   t_prev <- layer_prev$t
  #   
  #   degs_t_prev <- getDEGs(DEGs, "TM", t_prev)
  #   degs_t_prev <- degs_t_prev[degs_t_prev %in% V(gTfNetwork)$name]
  #   
  #   # Here we slipt them in two batches.
  #   batch1 <- intersect(degs_t_prev, degs_t) # code computed using t
  #   batch2 <- setdiff(degs_t_prev, batch1) # code computed using t-1
  #   
  #   # Take their 101, 111 and 011, and add them to the top layer
  #   batch1_codes <- sapply(batch1, getTemporalDECode, t)
  #   batch2_codes <- sapply(batch2, getTemporalDECode, t_prev)
  #   res$topLayer$`101` <- union(names(batch1_codes[batch1_codes == "101"]), res$topLayer$`101`)
  #   res$topLayer$`111` <- union(names(batch1_codes[batch1_codes == "111"]), res$topLayer$`111`)
  #   res$topLayer$`011` <- union(names(batch1_codes[batch1_codes == "011"]), res$topLayer$`011`)
  #   res$topLayer$`101` <- union(names(batch2_codes[batch2_codes == "101"]), res$topLayer$`101`)
  #   res$topLayer$`111` <- union(names(batch2_codes[batch2_codes == "111"]), res$topLayer$`111`)
  #   res$topLayer$`011` <- union(names(batch2_codes[batch2_codes == "011"]), res$topLayer$`011`)
  # }
  

  # Get all the degs active at time t-1 (version 2)
  t_prev <- NULL
  if (length(cascade) > 1) {
    layer_prev <- cascade[[length(cascade)]]
    t_prev <- layer_prev$t

    degs_t_prev <- getDEGs(DEGs, "TM", t_prev)
    degs_t_prev <- degs_t_prev[degs_t_prev %in% V(gTfNetwork)$name]

    # Takes the genes that have not already been attributed to one of
    # the three groups, i.e. remove the genes that are active at t
    # and are in one of the three groups buikt at t.
    alreadyExplained <- Reduce(union, c(
      res$topLayer$`101`,
      res$topLayer$`111`,
      res$topLayer$`011`
    ))
    degs_t_prev <- degs_t_prev[!degs_t_prev %in% alreadyExplained]
    # NOTE: we sill have genes that are active in t-1 and t.

    # Compute their code at t-1 and takes the genes that have the code 101, 111 or 011.
    degs_t_prev_codes <- sapply(degs_t_prev, getTemporalDECode, t_prev)

    # COMMENTED AS REQUESTED BY JENNY
    # 2017-03-13: Jenny wants a run of this script with the following lines uncommented (original version of the cascade)
    # ===== START =====
    res$topLayer$`101` <- union(names(degs_t_prev_codes[degs_t_prev_codes == "101"]), res$topLayer$`101`)
    res$topLayer$`111` <- union(names(degs_t_prev_codes[degs_t_prev_codes == "111"]), res$topLayer$`111`)
    res$topLayer$`011` <- union(names(degs_t_prev_codes[degs_t_prev_codes == "011"]), res$topLayer$`011`)
    # ===== END =====
  }

  
  # # To complete the top layer, we take all the active genes at t-1 that 
  # # had the code 001.
  # if (!is.null(t_prev)) {
  #   # degs_t_prev already set
  #   degs_t_prev_codes <- sapply(degs_t_prev, getTemporalDECode, t_prev)
  #   res$topLayer$past001 <- names(degs_t_prev_codes[degs_t_prev_codes == "001"])
  #   
  #   # Find all the past unexplained and substract them from the past 001.
  #   if (removeUnexplainedFromPast001) {
  #     pastUnexplained <- character(0)
  #     lastLayer <- NULL
  #     for (layer in cascade) {
  #       pastUnexplained <- union(layer$bottomLayer$groups$empty, pastUnexplained)
  #     }
  #     res$topLayer$past001 <- setdiff(res$topLayer$past001, pastUnexplained)
  #   }
  # }
  
  # To complete the top layer, we take all the active genes at t-1 that 
  # had the code 001 and *have been explained at t-1*.
  if (!is.null(t_prev)) {
    layer_t_prev <- cascade[[length(cascade)]]
    res$topLayer$past001 <- Reduce(union, c(
      layer_t_prev$bottomLayer$groups$`111`$degs,
      layer_t_prev$bottomLayer$groups$`101-011`$degs,
      layer_t_prev$bottomLayer$groups$`101-111-011`$degs,
      layer_t_prev$bottomLayer$groups$`001`$new$degs,
      layer_t_prev$bottomLayer$groups$`001`$past$degs
    ))
  }
  
  
  # Now we are ready to explain the *new* synergistic genes (001)
  new_degs_t <- newDegs[[paste0("TM_", t)]]
  new_degs_t <- new_degs_t[new_degs_t %in% V(gTfNetwork)$name]
  new_degs_t_codes <- sapply(new_degs_t, getTemporalDECode, t)
  
  res$bottomLayer$`001` <- names(new_degs_t_codes[new_degs_t_codes == "001"])
  toExplain <- res$bottomLayer$`001`
  
  # Rev 2e: instead of explaining the new 001 at t, we take *all* the 001 that are active at t.
  res$bottomLayer$`001` <- names(degs_t_codes[degs_t_codes == "001"])
  toExplain <- res$bottomLayer$`001`
  
  # Among the 001, we want to find the ones that can be explained by at least
  # one of the three motifs 111, 101-111-011 and 101-011.
  # First, we generate an igraph network with the degs of interest
  degs <- Reduce(union, c(
    res$topLayer$`101`,
    res$topLayer$`111`,
    res$topLayer$`011`,
    toExplain
  ))
  
  g <- delete.vertices(gTfNetwork, setdiff(V(gTfNetwork)$name, degs))
  
  # Finds the neighbors for each 001
  degs001_neighbors <- lapply(toExplain, getNeighbors, 1, g)
  names(degs001_neighbors) <- toExplain
  
  explained <- character(0)
  for (deg in toExplain) {
    neighbors <- degs001_neighbors[[deg]]
    
    hasNeighborIn101 <- (length(intersect(neighbors, res$topLayer$`101`)) > 0)
    hasNeighborIn111 <- (length(intersect(neighbors, res$topLayer$`111`)) > 0)
    hasNeighborIn011 <- (length(intersect(neighbors, res$topLayer$`011`)) > 0)
    
    if (hasNeighborIn111 && (!hasNeighborIn101 || !hasNeighborIn011)) { # S1
      res$bottomLayer$groups$`111`$degs <- union(deg, res$bottomLayer$groups$`111`$degs)
      # TODO: add edges
    } else if (hasNeighborIn101 && hasNeighborIn111 && hasNeighborIn011) { # S2
      res$bottomLayer$groups$`101-111-011`$degs <- union(deg, res$bottomLayer$groups$`101-111-011`$degs)
      # TODO: add edges
    } else if (hasNeighborIn101 && hasNeighborIn011 && !hasNeighborIn111) { # S3
      res$bottomLayer$groups$`101-011`$degs <- union(deg, res$bottomLayer$groups$`101-011`$degs)
      # TODO: add edges
    } else {
      # this 001 could still be explained with 2nd generation new 001
      # or past 001
      next
    }
    
    explained <- union(deg, explained)
  }
  # Update the list of genes to explain
  toExplain <- setdiff(toExplain, explained)
  
  
  # We will now try to explain the remaining genes using the three
  # sets of 001 that we have just computed (S1, S2, S3)
  # We start by making a new igraph that contains S1, S2, S3, and
  # the remaining 001.
  degs <- Reduce(union, c(
    res$bottomLayer$groups$`111`$degs,
    res$bottomLayer$groups$`101-111-011`$degs,
    res$bottomLayer$groups$`101-011`$degs,
    toExplain
  ))
  
  g <- delete.vertices(gTfNetwork, setdiff(V(gTfNetwork)$name, degs))
  
  # Finds the neighbors for the remaining 001
  degs001_neighbors <- lapply(toExplain, getNeighbors, 1, g)
  names(degs001_neighbors) <- toExplain
  
  explained <- character(0)
  for (deg in toExplain) {
    neighbors <- degs001_neighbors[[deg]]
    
    hasNeighborIn001_111 <- (length(intersect(neighbors, res$bottomLayer$groups$`111`$degs)) > 0)
    hasNeighborIn001_101_111_011 <- (length(intersect(neighbors, res$bottomLayer$groups$`101-111-011`$degs)) > 0)
    hasNeighborIn001_101_011 <- (length(intersect(neighbors, res$bottomLayer$groups$`101-011`$degs)) > 0)
    
    if (hasNeighborIn001_111 || hasNeighborIn001_101_111_011 || hasNeighborIn001_101_011) {
      res$bottomLayer$groups$`001`$new$degs <- union(deg, res$bottomLayer$groups$`001`$new$degs)
      # TODO: add edges
    } else {
      next
    }
    
    explained <- union(deg, explained)
  }
  # Update the list of genes to explain  
  toExplain <- setdiff(toExplain, explained)
  
  
  # We have the possibility to do a recursive search to explain
  # the remaining current 001 using current 001 that have already
  # been explained
  # DO NOT CLEAR THE EXPLAINED SET FROM THE PREVIOUS STEP
  lastStage <- explained
  if (recursiveSearch) {
    count <- 0
    repeat {
      if (length(lastStage) == 0 || length(toExplain) == 0) {
        break
      }
      count <- count + 1
      # cat("Counter time ", t, ": ", count, "\n", sep = "")
      explained <- character(0)
      for (deg in toExplain) {
        toKeep <- union(deg, lastStage)
        g2 <- delete.vertices(gTfNetwork, setdiff(V(gTfNetwork)$name, toKeep))
        neigh <- getNeighbors(deg, 1, g2)
        if (length(neigh) > 0) {
          res$bottomLayer$groups$`001`$new$degs <- union(deg, res$bottomLayer$groups$`001`$new$degs)
          # TODO: add edges
          
          explained <- union(deg, explained)
        }
      }
      toExplain <- setdiff(toExplain, explained)
      lastStage <- explained
    }
  }
  
  
  # We have tried to explain as many as possible new 001 using the motifs
  # 111, 101-111-011, 101-011 and new 001 already explained.
  # We will now try to explain the remaining new 001 using the past 001.
  if (!is.null(t_prev)) {
    degs <- Reduce(union, c(
      res$topLayer$past001,
      toExplain
    ))
    
    g <- delete.vertices(gTfNetwork, setdiff(V(gTfNetwork)$name, degs))
    
    # Finds the neighbors for the remaining 001
    degs001_neighbors <- lapply(toExplain, getNeighbors, 1, g)
    names(degs001_neighbors) <- toExplain
    
    # toExplain <- res$bottomLayer$`001`
    explained <- character(0)
    for (deg in toExplain) {
      neighbors <- degs001_neighbors[[deg]]
      
      hasNeighborInPast001 <- (length(intersect(neighbors, res$topLayer$past001)) > 0)
      
      if (hasNeighborInPast001) {
        res$bottomLayer$groups$`001`$past$degs <- union(deg, res$bottomLayer$groups$`001`$past$degs)
        # TODO: add edges
      } else {
        next
      }
      
      explained <- union(deg, explained)
    }
    # Update the list of genes to explain
    toExplain <- setdiff(toExplain, explained)
    
    
    # We have the possibility to do a recursive search to explain
    # the remaining 001 using the 001 previously explained using
    # the past 001.
    lastStage <- explained
    if (recursiveSearch) {
      count <- 0
      repeat {
        if (length(lastStage) == 0 || length(toExplain) == 0) {
          break
        }
        count <- count + 1
        # cat("Counter time ", t, ": ", count, "\n", sep = "")
        explained <- character(0)
        for (deg in toExplain) {
          toKeep <- union(deg, explained)
          g2 <- delete.vertices(gTfNetwork, setdiff(V(gTfNetwork)$name, toKeep))
          neigh <- getNeighbors(deg, 1, g2)
          if (length(neigh) > 0) {
            res$bottomLayer$groups$`001`$past$degs <- union(deg, res$bottomLayer$groups$`001`$past$degs)
            # we don't keep track of the edges
            
            explained <- union(deg, explained)
          }
        }
        toExplain <- setdiff(toExplain, explained)
        lastStage <- explained
      }
    }
  }
  
  
  # The remaining new 001 are considered "unexplained"
  res$bottomLayer$groups$empty <- toExplain
  
  
  # save
  cascade[[length(cascade)+1]] <- res
}


# layer <- cascade[[3]]
# layer$bottomLayer$groups$`001`$past$degs








# Prints some numbers about a layer of the cascade
printCascadeLayer <- function(cascadeLayer) {
  
  topLayer <- cascadeLayer$topLayer
  cat("====================\n")
  cat("Time: ", cascadeLayer$t, "\n", sep = "")
  cat("Num new degs: ", length(cascadeLayer$degs), "\n", sep = "")
  cat("\n")
  
  cat("Top layer:\n")
  cat("Number of 101 degs: ", length(topLayer$`101`), ": ", paste0(sort(topLayer$`101`), collapse = ", "), "\n", sep = "")
  cat("Number of 111 degs: ", length(topLayer$`111`), ": ", paste0(sort(topLayer$`111`), collapse = ", "), "\n", sep = "")
  cat("Number of 011 degs: ", length(topLayer$`011`), ": ", paste0(sort(topLayer$`011`), collapse = ", "), "\n", sep = "")
  cat("Number of past 001: ", length(topLayer$past001), "\n", sep = "")
  cat("\n")
  
  bottomLayer <- cascadeLayer$bottomLayer
  cat("Bottom layer:\n")
  cat("Number of 001 <- 111 degs: ", length(bottomLayer$groups$`111`$degs), ": ", 
      paste0(sort(bottomLayer$groups$`111`$degs), collapse = ", "), "\n", sep = "")
  # cat("Number of 001 <- 111 edges: ", nrow(bottomLayer$groups$`111`$undirEdges), "\n", sep = "")
  cat("\n")
  
  cat("Number of 001 <- (101 & 111 & 011) degs: ", length(bottomLayer$groups$`101-111-011`$degs), ": ",
      paste0(sort(bottomLayer$groups$`101-111-011`$degs), collapse = ", "), "\n", sep = "")
  # cat("Number of 001 <- 101 edges: ", nrow(bottomLayer$groups$`101-111-011`$undirEdges$source_101), "\n", sep = "")
  # cat("Number of 001 <- 111 edges: ", nrow(bottomLayer$groups$`101-111-011`$undirEdges$source_111), "\n", sep = "")
  # cat("Number of 001 <- 011 edges: ", nrow(bottomLayer$groups$`101-111-011`$undirEdges$source_011), "\n", sep = "")
  cat("\n")
  
  cat("Number of 001 <- (101 & 011) degs: ", length(bottomLayer$groups$`101-011`$degs), ": ",
      paste0(sort(bottomLayer$groups$`101-011`$degs), collapse = ", "), "\n", sep = "")
  # cat("Number of 001 <- 101 edges: ", nrow(bottomLayer$groups$`101-011`$undirEdges$source_101), "\n", sep = "")
  # cat("Number of 001 <- 011 edges: ", nrow(bottomLayer$groups$`101-011`$undirEdges$source_011), "\n", sep = "")
  cat("\n")
  
  # 001 that can be explained by new explained 001
  cat("Number of 001 <- new explained 001: ", length(bottomLayer$groups$`001`$new$degs), ": ",
      paste0(sort(bottomLayer$groups$`001`$new$degs), collapse = ", "), "\n", sep = "")
  # cat("Number of edges comming from 001 (111): ", nrow( bottomLayer$groups$`001`$new$`111`$undirEdges), "\n", sep = "")
  # cat("Number of edges comming from 001 (101-111-011): ", nrow( bottomLayer$groups$`001`$new$`101-111-011`$undirEdges), "\n", sep = "")
  # cat("Number of edges comming from 001 (101-011): ", nrow( bottomLayer$groups$`001`$new$`101-011`$undirEdges), "\n", sep = "")
  cat("\n")
  
  # 001 than can only be explained by past 001
  cat("Number of 001 <- past 001: ", length(bottomLayer$groups$`001`$past$degs), ": ",
      paste0(sort(bottomLayer$groups$`001`$past$degs), collapse = ", "), "\n", sep = "")
  # cat("Number of edges comming from past 001: ", nrow(bottomLayer$groups$`001`$past$undirEdges), "\n", sep = "")
  cat("\n")
  
  # 001 that cannot be explained
  cat("Number of 001 that cannot be explained: ", length(bottomLayer$groups$empty), ": ",
      paste0(sort(bottomLayer$groups$empty), collapse = ", "), "\n", sep = "")
}


layerIndex <- 6
layer <- cascade[[layerIndex]]
printCascadeLayer(cascade[[layerIndex]])

# sink(paste0("output/cascade_layer_", layerIndex, "_20160914.txt"))
# sink(paste0("output/cascade_layer_", layerIndex, "_20170128.txt"))
sink(paste0("output/cascade_layer_", layerIndex, "_20170313.txt"))
cat("top 101:", paste(layer$topLayer$`101`), "\n")
cat("top 111:", paste(layer$topLayer$`111`), "\n")
cat("top 011:", paste(layer$topLayer$`011`), "\n")
cat("past 001:", paste(layer$topLayer$past001), "\n")
cat("mid 111:", paste(layer$bottomLayer$groups$`111`$degs), "\n")
cat("mid 101-111-011:", paste(layer$bottomLayer$groups$`101-111-011`$degs), "\n")
cat("mid 101-011", paste(layer$bottomLayer$groups$`101-011`$degs), "\n")
cat("bot 001", paste(layer$bottomLayer$groups$`001`$new$degs), "\n")
cat("past 001 -> 001", paste(layer$bottomLayer$groups$`001`$past$degs), "\n")
cat("unexplained", paste(layer$bottomLayer$groups$empty))
sink()


fileConn <- file(paste0("output/cascade_layer_", layerIndex, "_20160914.txt"))
writeLines(c(paste0(layer$topLayer$`111`, sep = " "),
             "World"), fileConn)
close(fileConn)



# Looking for the edges (change index of layer to generate)
layer <- cascade[[1]]

# 111 -> 001
source <- layer$topLayer$`111`
target <- layer$bottomLayer$groups$`111`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])

source <- layer$topLayer$`101`
target <- layer$bottomLayer$groups$`101-111-011`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])
source <- layer$topLayer$`111`
target <- layer$bottomLayer$groups$`101-111-011`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])
source <- layer$topLayer$`011`
target <- layer$bottomLayer$groups$`101-111-011`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])

source <- layer$topLayer$`101`
target <- layer$bottomLayer$groups$`101-011`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])
source <- layer$topLayer$`011`
target <- layer$bottomLayer$groups$`101-011`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])

#
source <- layer$bottomLayer$groups$`111`$degs
target <- layer$bottomLayer$groups$`001`$new$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])
source <- layer$bottomLayer$groups$`101-111-011`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])
source <- layer$bottomLayer$groups$`101-011`$degs
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])

#
length(source <- layer$topLayer$past001)
length(target <- layer$bottomLayer$groups$`001`$past$degs)
nrow(tfNetwork[tfNetwork$Regulator %in% source & tfNetwork$Target %in% target,])

sum(source %in% target)

# # ====================
# # Converts Eren's viper file
# 
# library(compositions)
# 
# filename <- "data/viper/motif_analysis.csv"
# viper <- read.csv(filename, header = TRUE)
# rownames(viper) <- viper[,1]
# viper[,1] <- NULL
# 
# viperFinal <- matrix(0, nrow=nrow(viper), ncol=ncol(DEGs$data))
# rownames(viperFinal) <- rownames(viper)
# colnames(viperFinal) <- colnames(DEGs$data)
# 
# for (i in 1:nrow(viperFinal)) {
#   gene <- rownames(viperFinal)[i]
#   for (t in timePoints) {
#     if (viper[gene, paste0("Viper_sumarry_", t)] > 7)
#       stop("Larger than expected max value = 7")
#     code <- binary(viper[gene, paste0("Viper_sumarry_", t)]) # "sumarry" in the file instead of "summary"
#     code <- paste0(paste0(rep("0", 3-nchar(code)), collapse = ""), code) # pad with 0
#     # cat(code, ": T=", strtoi(substr(code, 1, 1)), ", M=", strtoi(substr(code, 2, 2)), ", TM=", strtoi(substr(code, 3, 3)), "\n", sep = "")
#     
#     # keep only the last three digits
#     code <- substr(code, nchar(code)-2, nchar(code))
#     
#     viperFinal[gene, paste0("T_", t)] <- - strtoi(substr(code, 1, 1))
#     viperFinal[gene, paste0("M_", t)] <- - strtoi(substr(code, 2, 2))
#     viperFinal[gene, paste0("TM_", t)] <- - strtoi(substr(code, 3, 3))
#   }
# }
# # saveRDS(viperFinal, "data/viper/viper_activity_data.rds")
# 
# # VIPERs <- list(
# #   name = "Viper",
# #   data = readRDS("data/viper/viper_activity_data.rds"),
# #   threshold = -1e-20
# # )
