library(tidyverse)
library(TopDom)
library(strawr)
library(hicrep)
set.seed(2022)

hic2matrix <- function(file_path, chrom, output_chr = chrom, resol, method = "NONE") {
  hic_mat <-
    hicrep::hic2mat(
      file = file_path,
      chromosome1 = chrom,
      chromosome2 = chrom,
      resol = resol,
      method = method
    )

  # Build header columns
  col_start <- as.integer(seq(0, (nrow(hic_mat) - 1) * resol, by = resol))
  col_end <- col_start + as.integer(resol)

  dt <- bind_cols(
    tibble(chrom = output_chr, start = col_start, end = col_end),
    as_tibble(hic_mat)
  )

  output_path <- tempfile(fileext = ".txt")
  on.exit(unlink(output_path), add = TRUE)

  write_delim(dt, output_path, na = "0", delim = " ", col_names = FALSE)
  TopDom::readHiC(output_path)
}



args <- commandArgs(trailingOnly = TRUE)
hic_file <- args[1]
output_domain <- args[2]

logging::loginfo(str_interp("${hic_file} : ${output_domain}"))

topdom_results <- c(as.character(1:22), "X") %>%
  map(function(chrom) {
    logging::loginfo(str_interp("Processing ${hic_file}: chr${chrom} ..."))
    hic_matrix <- hic2matrix(file_path = hic_file,
      chrom = chrom, resol = 50e3L, method = "KR")
    TopDom::TopDom(hic_matrix, window.size = 5, debug = TRUE)
  })


map_dfr(topdom_results, function(item) {
  item$bed %>% as_tibble() %>% filter(name == "domain") %>%
  select(-name) %>%
  mutate(chrom = as.character(chrom)) %>%
  mutate(chrom2 = chrom, chromStart2 = chromStart, chromEnd2 = chromEnd)
}) %>%
  write_tsv(output_domain)


