# ENV
rm(list = ls())
gc()

library(tidyverse)
library(stringr)
library(rjson)

# wk <- "D:\\WeChatAfficialAccount\\KEGG-database\\第二版\\kegg_compound"
wk <- "/Users/xiaofei/Desktop/KEGG"
setwd(wk)


# Parsed Json
cpd <- fromJSON(file = "br08901.json", simplify = TRUE)
level <- cpd$name
json <- cpd$children

rlt <- data.frame()

pbapply::pblapply(1:length(json), FUN = function(a){
  level1 <- json[[a]]$name
  json1 <- json[[a]]$children
  
  for (b in 1:length(json1)){
    level2 <- json1[[b]]$name
    json2 <- json1[[b]]$children
    
    if (!is.null(json2)){
      for (c in 1:length(json2)){
        level3 <-json2[[c]]$name
        term <- data.frame(
          "BRID" = "br08901",
          "L1" = level1,
          "L2" = level2,
          "mapID" = paste0("map",str_extract(level3, pattern = "([0-9]{5})(\\s+?)(.*)", group = 1)),
          "mapName" = str_extract(level3, pattern = "([0-9]{5})(\\s+?)(.*)", group = 3)
        )
        rlt <<- rbind(rlt, term)
      }
    }else{
      message("++++++++++++++++++++++")
      message(level1)
      level3 <- NA
      term <- data.frame(
        "BRID" = "br08901",
        "L1" = level1,
        "L2" = NA,
        "L3" = NA
      )
      rlt <<- rbind(rlt, term)
    }
  }
})

# pathway to compound
library(magrittr)
map2cpd <- read.table("pathway2compound.txt", header = FALSE, sep = "\t", quote = "") %>%
  `colnames<-`(c("mapID", "CPD"))
map2cpd %<>% mutate(
  mapID = map_chr(mapID, \(x) strsplit(x, ":", fixed = TRUE)[[1]][2]),
  CPD = map_chr(CPD, \(x) strsplit(x, ":", fixed = TRUE)[[1]][2])
)

# compound information
cpd2info <- read.table("compound.txt", header = FALSE, sep = "\t", quote = "") %>%
  `colnames<-`(c("CPD", "CPDName"))



cls <- left_join(x = rlt, y = map2cpd, by = "mapID")

cls2 <- left_join(x = cls, y = cpd2info, by = "CPD")












