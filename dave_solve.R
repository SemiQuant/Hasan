# RULES
# (1) All information for each element is shown in one sequence or another (ie there is no chance of a letter appearing 
# in a position in a sequence if it hasn't appeared in that position in a sequence where that element is revealed)
# (2) Only sequence connections that are revealed are assumed to be possible (such as between element 4 and 5), 
# unless there is no sequence connection at all, in which case any of the revealed elements are possible connections (ie between element 6 and 7)

# rule 2 not implemented in the below




# possible_characters <- c("A", "G", "C", "T", "–")
# NA = no information
dat <- tibble::tribble(
  ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`,
  1001L,  "A",  "A",  "A",  "G",  NA,  NA,  NA,  NA,  NA,   NA,
  9020L,  "A",  "A",  "A",  "C",  NA,  NA,  NA,  NA,  NA,   NA,
  506L,  NA,  NA,  NA,  "G",  "A",  "A",  NA,  NA,  NA,   NA,
  4035L,  NA,  NA,  NA,  "C",  "T",  "A",  NA,  NA,  NA,   NA,
  630L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "–",   "G",
  4025L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "A",   "G"
)
# 
# dat <- tibble::tribble(
#   ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`,
#      10L,  "X",  "A",  "A",  "G",  "A",  "A",  "N",  "N",  "N",   "N",
#      90L,  "Y",  "A",  "A",  "C",  "T",  "A",  "N",  "N",  "N",   "N",
#      10L,  "N",  "N",  "N",  "N",  "N",  "N",  "A",  "C",  "X",   "G",
#      90L,  "N",  "N",  "N",  "N",  "N",  "N",  "A",  "C",  "Y",   "G"
#   )
# dat[dat == "N"] <- NA
# 






require(tidyverse)
rl <- as.relistable(sapply(dat[,-1], function(x){unique(na.omit(x))}))
r <- expand.grid(data.frame(rl, stringsAsFactors = F), stringsAsFactors = F) %>% 
  distinct()
frequency_assignment <- matrix(0, nrow = nrow(r), ncol = nrow(dat))
for(x in 1:nrow(dat)){
  valid <- apply(r, 1, function(z){
    all(z == dat[x, -1], na.rm = T)
  })
  frequency_assignment[valid , x] = dat$Count[x]/sum(valid)
}

sequences <- data.frame(frequency = rowSums(frequency_assignment),
                       probability = rowSums(frequency_assignment)/sum(frequency_assignment),
                       r)
sequences <- sequences[sequences$frequency>0,]
sequences
# 
# need=2
# for (c in 3:ncol(sequences)){
#   if (n_distinct(sequences[,c]) == need){
#     sequences <- sequences %>% 
#       group_by(sequences[,c]) %>% 
#       filter(probability == min(probability)) %>% 
#       # slice(1) %>% # takes the first occurrence if there is a tie
#       ungroup() %>% 
#       mutate(frequency = frequency/sum(frequency)) %>% 
#       select(-ncol(.))
#       # select(cost, everything())
#     break
#   }
# }
# sequences
# 
# 
