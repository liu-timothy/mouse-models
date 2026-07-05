# install.packages("pamr")
library(pamr)

# ----------------------------
# User input:
# rows = genes, columns = samples
# first column should be GeneSymbol
# expression values should be VST/rlog/log2(TPM+1), not DESeq2 results
# ----------------------------
expr_file <- "rnaseq_expression_matrix.tsv"
out_file  <- "BASE47_predictions.tsv"

# ----------------------------
# Download official RNA-seq BASE47 predictor
# ----------------------------
predictor_url <- "https://raw.githubusercontent.com/kimlabunc/Kardos_Rose_NanoString/main/BASE47_RNASEQ_PREDICTOR_DO_NOT_EDIT.txt"
predictor_file <- "BASE47_RNASEQ_PREDICTOR_DO_NOT_EDIT.txt"

if (!file.exists(predictor_file)) {
  download.file(predictor_url, predictor_file)
}

# ----------------------------
# Load BASE47 predictor
# ----------------------------
x <- read.table(
  predictor_file,
  sep = "\t",
  row.names = 1,
  header = TRUE,
  check.names = FALSE
)

classes <- as.vector(t(x[1, ]))

xn <- apply(x[-1, ], 2, function(z) as.numeric(as.vector(z)))
rownames(xn) <- rownames(x)[-1]

base47_genes <- rownames(xn)

# ----------------------------
# Load your expression matrix
# ----------------------------
expr <- read.table(
  expr_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(expr) <- expr$GeneSymbol
expr$GeneSymbol <- NULL

# Keep protein-coding only if you have already filtered elsewhere.
# Not required here because BASE47 gene symbols define the needed genes.

# Collapse duplicate gene symbols by mean
expr <- aggregate(. ~ rownames(expr), data = expr, FUN = mean)
rownames(expr) <- expr[, 1]
expr <- expr[, -1]

# Make numeric matrix
expr <- as.matrix(expr)
mode(expr) <- "numeric"

# ----------------------------
# Check required genes
# ----------------------------
missing_genes <- setdiff(base47_genes, rownames(expr))

if (length(missing_genes) > 0) {
  stop(
    paste0(
      "Missing BASE47 genes: ",
      paste(missing_genes, collapse = ", ")
    )
  )
}

# Put genes in same order as predictor
test <- expr[base47_genes, ]

# ----------------------------
# Train and predict, following official code
# ----------------------------
trainSet <- list(
  x = scale(xn),
  y = classes,
  geneid = rownames(xn),
  genenames = rownames(xn)
)

mytrain <- pamr.train(trainSet)

pred.class <- pamr.predict(mytrain, scale(test), threshold = 0)
pred.prob  <- pamr.predict(mytrain, scale(test), threshold = 0, type = "posterior")

results <- data.frame(
  Sample = colnames(test),
  Basal_score = pred.prob[, 1],
  Luminal_score = pred.prob[, 2],
  BASE47_subtype = ifelse(pred.class == 1, "Basal", "Luminal")
)

write.table(
  results,
  out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(results)