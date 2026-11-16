# Export the Mayo Clinic primary biliary cholangitis (PBC) dataset from the

# Run from the project root:  Rscript code/00_export_data.R
suppressPackageStartupMessages(library(survival))
data(pbc, package = "survival")
write.csv(pbc, "C:/Users/User/Desktop/Gawesha/pbc-liver-survival-r/pbc-liver-survival-r/data/pbc.csv", row.names = FALSE, na = "")
cat(sprintf("Wrote data/pbc.csv: %d rows x %d columns\n", nrow(pbc), ncol(pbc)))

