# 
library(plyr)
test_df <- read.delim("clipboard", header = T, check.names = F)
y <- count(test_df, 'annotation')
