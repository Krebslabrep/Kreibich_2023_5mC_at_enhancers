######### Functions for Single molecule analysis #########

# Coverage matrix ----------------------------------------------------------------------------------------------------
## input is a CGme matrix list with 
## 1. level = sample - should be named !!!
## 2. level = bin
## 3. level = fractions (0|1)
## 4. level = methylation of each single molecule (single methylation for CpG, mean methylation within bin for GpC)

  make.coverage.tibble <- function(me_list){
    count_bin_l <- mclapply(seq(me_list), function(sample){
      # sample = 1
      bins_l <- mclapply(seq(me_list[[sample]]), function(bin){
        # bin = 12
        fractions_l <- lapply(seq(me_list[[sample]][[bin]]), function(fraction){
          # fraction = 1
          as.tibble(me_list[[sample]][[bin]][[fraction]][!is.na(me_list[[sample]][[bin]][[fraction]])]) %>%
            count(value)
        })
        names(fractions_l) <- names(me_list[[sample]][[bin]])
        bind_rows(fractions_l, .id = "fraction")
      }, mc.cores = 2)
      names(bins_l) <- names(me_list[[sample]])
      bind_rows(bins_l, .id = "bin")
    }, mc.cores = length(SAMPLENAMES))
    names(count_bin_l) <- names(me_list)
    count_bin_tbl <- bind_rows(count_bin_l, .id = "sample")
    return(count_bin_tbl)
  }
  
  
# CMH matrix ----------------------------------------------------------------------------------------------------
## input is a CGme matrix list with 
## 1. level = sample - should be named !!!
## 2. level = bin
## 3. level = fraction (0|1)
## 4. level = methylation of each single moelcule (single methylation for CpG, mean methylation within bin for GpC)
  
  filter.for.coverage <- function(CGme_list_filtered, cO, SAMPLENAMES, TF = FALSE){
      
    ## Filter for coverage and set remove bins with only one fraction - for every sample
      CGme_list2 <- mclapply(seq_along(CGme_list_filtered), function(sample){
        # sample = 1
        CGme_sample <- CGme_list_filtered[[sample]]
        CGme_sample2 <- mclapply(seq_along(CGme_sample), function(bin){
          # bin = 107
          if(TF == TRUE){fraction_names <- names(CGme_sample[[bin]])}
          CGme_sample3 <- lapply(seq_along(CGme_sample[[bin]]), function(fraction){
            # fraction = 1
            CGme_sample[[bin]][[fraction]] <- CGme_sample[[bin]][[fraction]][!is.na(CGme_sample[[bin]][[fraction]])]
            if(!length(CGme_sample[[bin]][[fraction]]) == 0){
              return(CGme_sample[[bin]][[fraction]])
            } else {
              CGme_sample[[bin]][[fraction]] <- NULL
            }
          })
          if(TF == TRUE){
            names(CGme_sample3) <- fraction_names
            if(length(CGme_sample3) > 2){
              if(length(unlist(CGme_sample3)) < cO){
                CGme_sample3 <- NULL
              } else {
                CGme_sample3
              }
            } else {
              CGme_sample3 <- NULL
            }
          } else {
            if(length(CGme_sample3) == 2){
              if(sum(length(CGme_sample3[[1]]), length(CGme_sample3[[2]])) < cO){
                CGme_sample3 <- NULL
              } else {
                CGme_sample3
              }
            } else {
              CGme_sample3 <- NULL
            }
          }
        }, mc.cores = 3)
        names(CGme_sample2) <- names(CGme_sample)
        return(CGme_sample2)
      }, mc.cores = 10)
      names(CGme_list2) <- SAMPLENAMES
      # View(CGme_list2)
      message("Filter for coverage complete")
      
      return(CGme_list2)
  }   
  
  make.CMH.input.tibble <- function(CGme_list2, SAMPLENAMES, BINNAMES){
  
    ## Combine for each bin the fraction CG methylation objects of all samples
      CGme_bin_list <- lapply(seq_along(CGme_list2[[1]]), function(bin){
        # print(bin)
        # bin = 34210
        if(length(SAMPLENAMES) == 3){
          bin_sample_list <- list(CGme_list2[[1]][[bin]], CGme_list2[[2]][[bin]], CGme_list2[[3]][[bin]])
        } else if (length(SAMPLENAMES) == 2){
          bin_sample_list <- list(CGme_list2[[1]][[bin]], CGme_list2[[2]][[bin]])
        }
        names(bin_sample_list) <- SAMPLENAMES
        return(bin_sample_list)
      })
      names(CGme_bin_list) <- BINNAMES
      # View(CGme_bin_list)
      message("Combining methylation matrixes complete")
    
    ## Combine lists to tibble containing with columns: Sample, Fraction, Value (Methylation)
      bin_tbl_list <- mclapply(seq_along(CGme_bin_list), function(bin){
        # bin = 1
        # message(bin)
        tl2 <- lapply(seq_along(CGme_bin_list[[bin]]), function(sample){
          # sample = 2
          if(!is.null(CGme_bin_list[[bin]][[sample]])){
            plyr::compact(CGme_bin_list[[bin]][[sample]])
            tl1 <- lapply(seq_along(CGme_bin_list[[bin]][[sample]]), function(fraction){
              tl0 <- as_tibble(CGme_bin_list[[bin]][[sample]][[fraction]]) %>%
                drop_na()
              rbind(tl0,  tibble(value = c(0,1))) #add pseudo count of 1 for both values
            })
          } else { 
            tl1 <- NULL 
          }
          names(tl1) <- names(CGme_bin_list[[bin]][[sample]])
          bind_rows(tl1, .id = "fraction")
        })
        names(tl2) <- names(CGme_bin_list[[bin]])
        
        tbl <- bind_rows(tl2, .id = "sample")
        if(nrow(tbl) != 0) {
          tbl <-  tbl %>% filter(!is.na(value))
        } else {
          tbl <- NULL
        }
        tbl
      }, mc.cores = 10)
      names(bin_tbl_list) <- BINNAMES
      # View(bin_tbl_list)
      message("Making tibbles for all bins complete")
      
      return(bin_tbl_list)
  }
  
  make.CMH.input.tibble.TF <- function(CGme_list2, SAMPLENAMES){
    ## Convert to tibbles and filter for TFBS bin that contain bound fraction
    CGme_tbl_l <- lapply(seq_along(CGme_list2), function(sample){
      bin_tbl_l <- mclapply(seq_along(CGme_list2[[sample]]), function(bin){
        fraction_tbl_l <- lapply(seq_along(CGme_list2[[sample]][[bin]]), function(fraction){
          tbl1 <- as_tibble(CGme_list2[[sample]][[bin]][[fraction]])
          rbind(tbl1,  tibble(value = c(0,1))) #add pseudo count of 1 for both values
        })
        names(fraction_tbl_l) <- str_remove(names(CGme_list2[[sample]][[bin]]), "TFBS_.{1,}_")
        if(sum(grepl("101", names(fraction_tbl_l))) != 1){
          fraction_tbl_l <- NULL
        } else {
          bind_rows(fraction_tbl_l, .id = "fraction")
        }
      }, mc.cores = 10)
      names(bin_tbl_l) <- names(CGme_list2[[sample]])
      bind_rows(compact(bin_tbl_l), .id = "TFBS")
    })
    names(CGme_tbl_l) <- SAMPLENAMES
    message("Converting to tibbles complete")
    return(CGme_tbl_l)
  }
  
# CMH testing ----------------------------------------------------------------------------------------------------
## input is a list of tibbles with 
## 1. level = bin - should be named!! Is NULL if not covered in at least two replicates with a coverage cutoff cO
## 2. tibble with sample, fraction, value -->  methylated read counts per fraction per sample

  
  make.3Dcontig.table <- function(bin_tbl_list, BINNAMES){
    
    ## Compute 3dimensional contingency tables in an array format for each bin
      bin_arrays <- mclapply(seq_along(bin_tbl_list), function(bin){
        # message(bin)
        if(!is.null(bin_tbl_list[[bin]])){
          table(bin_tbl_list[[bin]]$fraction, bin_tbl_list[[bin]]$value, bin_tbl_list[[bin]]$sample)
        } else {
          NULL
        }
      }, mc.cores = 10)
      # View(bin_arrays)
      names(bin_arrays) <- BINNAMES
      message("Making 3D contingency arrays complete")
      return(bin_arrays)
  }

  make.3Dcontig.table.TF <- function(freq_tbl, BINNAMES){
    
    ## Compute 3dimensional contingency tables in an array format for each bin
      bin_arrays <- mclapply(BINNAMES, function(BINNAME){
        TFBS_tbl <- freq_tbl %>% 
          filter(TFBS == BINNAME,
                 state != "accessible") %>% 
          mutate(state = factor(state, levels = c("closed", "bound"))) %>% 
          mutate(value = replace(value, value == 0.5, 0))
        table(TFBS_tbl$state, TFBS_tbl$value, TFBS_tbl$sample)
      }, mc.cores = 10)
      message("Making 3D contingency arrays complete")
      return(bin_arrays)
  }
  
  
  perform.CMH.test <- function(bin_arrays, BINNAMES){

    ## Perform Cochran-Mantel-Haenszel test 
      CMH_list <- mclapply(seq_along(bin_arrays), function(bin){
        # message(bin)
        # bin = 1
        if(!is.null(bin_arrays[[bin]])){
          if(dim(bin_arrays[[bin]])[1] == 2 & dim(bin_arrays[[bin]])[2] == 2 & (dim(bin_arrays[[bin]])[3] >= minSample)){
            mantelhaen.test(bin_arrays[[bin]])
          }
        } 
      }, mc.cores = 10)
      # View(CMH_list)
      message("Performing Cochran-Mantel-Haenszel test complete")
    
    ## Extract pvalue and odds ratio from Cochran-Mantel-Haenszel test 
      CMH_pval <- lapply(seq_along(CMH_list), function(bin){
        # bin = 1
        c("pval" = CMH_list[[bin]]$p.val, CMH_list[[bin]]$estimate)#, CMH_list[[bin]]$statistic)
      })
      names(CMH_pval) <- BINNAMES
      # View(CMH_pval)
    
    ## Remove all empty bins
    CMH_pval <- plyr::compact(CMH_pval)
    BINNAMES_compact <- names(CMH_pval)
    
    ## Bind pvalue and odds ratio of all bins together into one tibble
    Final_tbl <- bind_rows(CMH_pval, .id = "bin") %>% 
      dplyr::rename(COR = `common odds ratio`) %>%
      mutate(bin = BINNAMES_compact)
 
    return(Final_tbl)
  }
  

  
  # Combine tibbles from all chromosomes to one big tibble -----------------------------------
  
  combine.chr.tibbles <- function(OUTDIR, OUTPUTNAME){

    files <- list.files(path = OUTDIR, pattern = paste0("^", OUTPUTNAME, '_chr\\d{1,2}'))
    if(length(files) == 19){
      CMH_tbl_l <- lapply(seq(19), function(CHR){
        # message(CHR)
        input_file <- paste0(OUTDIR, OUTPUTNAME, '_chr', CHR, '.rds')
        l <- readRDS(input_file)
        # file.remove(input_file)
        return(l)
      })
      if(class(CMH_tbl_l[[1]])[1] == "character"){
          CMH_tbl <- unlist(CMH_tbl_l)
        } else { 
          CMH_tbl <- bind_rows(CMH_tbl_l)
        }
      
      saveRDS(CMH_tbl, paste0(OUTDIR, OUTPUTNAME, '.rds'))
      message(paste("Saving", paste0(OUTDIR, OUTPUTNAME, '.rds'),"complete"))
      
    } else {
      x <- 19-length(files)
      message(paste(x, "files still missing of ", OUTPUTNAME))
    }
  }
  