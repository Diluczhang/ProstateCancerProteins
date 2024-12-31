# Function: Create qs.results and put it into a RDS file
# Input:    cluster (cluster name)
ProcessQSresults <- function (cluster) {
  print(crtwd)
  outwd <- paste(task.path, "QuSAGE_rds", label_gmt, sep = "/")
  if (!file.exists(outwd)) dir.create(outwd, recursive = T)
  file_qs_result <- paste(outwd, "/", "Cluster_", cluster, ".qs.rds", sep = "") #print(file_qs_result)
  if (!file.exists(file_qs_result)) {
    print(length(eset)); print(length(labelset));
    qs.results <<- qusage(eset, labelset, contrast, genesets)
    saveRDS(qs.results, file = file_qs_result)
  } else {
      print(file_qs_result)
      qs.results <<- readRDS(file_qs_result)
  }
}