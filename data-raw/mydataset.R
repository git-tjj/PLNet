## code to prepare `mydataset` dataset goes here

realdatamatrix<-readRDS(file = "data-raw/realdatamatrix.Rdata")
obs_mat<-readRDS(file = "data-raw/obs_mat.Rdata")
usethis::use_data(obs_mat, overwrite = FALSE)
usethis::use_data(realdatamatrix, overwrite = FALSE)
