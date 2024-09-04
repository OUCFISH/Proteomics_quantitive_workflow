rm(list = ls())
suppressWarnings(suppressMessages(library(openxlsx)))
sample.info <- read.xlsx(Sys.glob(paste0("data/","*分组*.xlsx")),startRow = 11)
sample.data <- sample.info[,c(1,2,3)]
# 定义填充函数
fill_na <- function(vec) {
  last_value <- NA
  for (i in 1:length(vec)) {
    if (!is.na(vec[i])) {
      last_value <- vec[i]
    } else {
      vec[i] <- last_value
    }
  }
  return(vec)
}
sample.data$X3 <- fill_na(sample.data$X3)
sample.data <- na.omit(sample.data)
colnames(sample.data) <- c("SampleID","SampleName","Group")
write.xlsx(x = sample.data,file = "./temp/group.xlsx")
## contrast
contrast <- na.omit(sample.info[c("对照组","处理组")])
contrast.name <- data.frame(contrast = paste0(contrast$对照组,"_vs_",contrast$处理组))
write.xlsx(x = contrast.name,file = "./temp/contrast.xlsx")

