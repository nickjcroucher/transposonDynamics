#!/usr/bin/env Rscript
# combine_prot_json.R
# Stacks every prot_NNNN.json page from the InterPro API into one data frame.
#   Shell : Rscript combine_prot_json.R [json_dir]   (default json_dir: prot_json)
#   R     : source("combine_prot_json.R")            # leaves `prot_df` in your session
# Needs : jsonlite (its rbind_pages() uses vctrs under the hood)

library(jsonlite)

args     <- commandArgs(trailingOnly = TRUE)
json_dir <- if (length(args) > 0) args[1] else "prot_json"

# 1. Page files, in numeric order
files <- list.files(json_dir, pattern = "^prot_[0-9]+\\.json$", full.names = TRUE)
if (length(files) == 0) stop("No prot_*.json files found in '", json_dir, "'")
files <- files[order(as.integer(gsub("[^0-9]", "", basename(files))))]

# 2. Read each page; flatten = TRUE turns nested objects into columns
#    (metadata.accession, metadata.name, ..., extra_fields.counters.proteins)
read_page <- function(f) {
  res <- fromJSON(f, flatten = TRUE)$results
  if (!is.data.frame(res) || nrow(res) == 0) return(NULL)
  res$page_file <- basename(f)          # lets you trace any row back to its page
  res
}
pages <- Filter(Negate(is.null), lapply(files, read_page))
if (length(pages) == 0) stop("The JSON files contain no results")

# 3. Stack all pages (a column missing from a page is filled with NA)
prot_df <- rbind_pages(pages)

# 4. Sanity checks against the total the API reports
api_count <- fromJSON(files[1])$count
if (is.null(api_count)) api_count <- NA_integer_
message(sprintf("%d files -> %d rows (API count: %s); %d duplicated accessions",
                length(files), nrow(prot_df), format(api_count),
                sum(duplicated(prot_df$metadata.accession))))
if (!is.na(api_count) && nrow(prot_df) != api_count) {
  warning("Row count differs from the API count: a page may be missing, or the ",
          "data changed during the download (e.g. a new InterPro release).")
}

# 5. Optional tidy-up: shorter column names, trimmed lineage strings
names(prot_df) <- sub("^metadata\\.", "", names(prot_df))
names(prot_df) <- sub("^extra_fields\\.counters\\.", "n_", names(prot_df))
if ("lineage" %in% names(prot_df)) prot_df$lineage <- trimws(prot_df$lineage)

str(prot_df)

# 6. Keep a copy on disk if you like (uncomment one line)
# saveRDS(prot_df, "prot_combined.rds")
write.csv(prot_df, "prot_combined.csv", row.names = FALSE)
