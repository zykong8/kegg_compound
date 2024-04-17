# env
rm(list = ls())
gc()

# direcroty
setwd("D:\\WeChatAfficialAccount\\KEGG-database")

library(tidyverse)
library(stringr)
library(rjson)

# All compounds
GetCompounds <- function(brs){
    # Save
    rlt <- data.frame()
    for (br in brs){
        # Download
        if (!file.exists(paste0("kegg_compound_",br,".json"))){
          message(paste0("Dwonload:", br, ", and parse JSON!"))
          src <- paste0("https://rest.kegg.jp/get/br:",br,"/json")
          download.file(url = src, destfile = paste0("kegg_compound_",br,".json"))
        }
        
        # Parsed Json
        cpd <- fromJSON(file = paste0("kegg_compound_",br,".json"), simplify = TRUE)
        level <- cpd$name
        json <- cpd$children
    
        # Extract compounds
        
        pbapply::pblapply(1:length(json), FUN = function(a){
            level1 <- json[[a]]$name
            json1 <- json[[a]]$children
            
            for (b in 1:length(json1)){
              level2 <- json1[[b]]$name
              json2 <- json1[[b]]$children
              
              if (!is.null(json2)){
                for (c in 1:length(json2)){
                  level3 <-json2[[c]]$name
                  json3 <- json2[[c]]$children
                  
                  if (!is.null(json3)){
                    for (d in 1:length(json3)){
                      level4 <- json3[[d]]$name
                      term <- data.frame(
                        "BRID" = br,
                        "L1" = level1,
                        "L2" = level2,
                        "L3" = level3,
                        "CID" = str_extract(level4, pattern = "(C[0-9]{5}) (.*)", group = 1),
                        "Compound" = str_extract(level4, pattern = "(C[0-9]{5}) (.*)", group = 2)
                      )
                      rlt <<- rbind(rlt, term)
                    }
                  }else{
                    level4 <- NA
                    term <- data.frame(
                      "BRID" = br,
                      "L1" = level1,
                      "L2" = level2,
                      "L3" = NA,
                      "CID" = str_extract(level3, pattern = "(C[0-9]{5}) (.*)", group = 1),
                      "Compound" = str_extract(level3, pattern = "(C[0-9]{5}) (.*)", group = 2)
                    )
                    rlt <<- rbind(rlt, term)
                  }
                }
              }else{
                level3 <- NA
                level4 <- NA
                term <- data.frame(
                  "BRID" = br,
                  "L1" = level1,
                  "L2" = NA,
                  "L3" = NA,
                  "CID" = str_extract(level2, pattern = "(C[0-9]{5}) (.*)", group = 1),
                  "Compound" = str_extract(level2, pattern = "(C[0-9]{5}) (.*)", group = 2)
                )
                rlt <<- rbind(rlt, term)
              }
            }
          }
        )
    }
    return(rlt)
}

brs <- c("br08001", "br08002", "br08003", "br08005",
         "br08006", "br08007", "br08009", "br08021")
cpds <- GetCompounds(brs = brs)
write.table(x = cpds, file = "kegg_compound_class.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)

# Create KEGG Compound GMT file
## Step 1. Download file
src2 <- "https://rest.kegg.jp/link/compound/pathway"
download.file(url = src2, destfile = "kegg_compound_pathway.txt")

src3 <- "https://rest.kegg.jp/list/pathway"
download.file(url = src3, destfile = "kegg_pathway_description.txt")

## Step 2. Extract information
map2cpd <- read.table("kegg_compound_pathway.txt", col.names = c("MapID", "CPD"))
map2cpd$CPD <- lapply(map2cpd$CPD, function(cid){
  str_remove(cid, "cpd:")
}) %>% unlist
map2cpd$MapID <- lapply(map2cpd$MapID, function(cid){
  str_remove(cid, "path:")
}) %>% unlist

map2cpd2tbl <- map2cpd %>% group_by(MapID) %>%
  reframe(
    "CID" = paste0(CPD, collapse = "\t")
  )

map2des <- read.table("kegg_pathway_description.txt", 
                      col.names = c("MapID", "Description"),
                      sep = "\t")
map2des$http <- paste0(c(src2, src3), collapse = ";")

## Step 3. Merge datasets
map2cid2des <- left_join(x = map2des, y = map2cpd2tbl, by = "MapID") %>%
  filter(!is.na(CID))

## Step 4. Save GMT file
write.table(x = map2cid2des[,2:4], file = "kegg_compound.gmt", 
            row.names = FALSE, col.names = FALSE,
            sep = "\t", quote = FALSE)