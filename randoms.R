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

for (c in 2:ncol(dat)){
  dat[,c] <- factor(as.character(unlist(dat[,c])), 
                    levels = c("A", "G", "C", "T", "–"))
}

require(mltools)
require(data.table)
one_hot(as.data.table(dat[-c(1)]))



dat
for (c in 2:ncol(dat)){
  dat[,c] <- paste0(unlist(dat[,c]), "_", colnames(dat)[c])
}
for (r in 1:nrow(dat)){
  dat[r,][grepl("NA_", dat[r,])] <- NA
}

test <- NULL
for (c in seq(ncol(dat)-1)){
  if (c>1){
    test.tmp <- dat[,c(1, c, c+1)]
    colnames(test.tmp) <- c("weight", "from", "to")
    test <- rbind(test, test.tmp)
  }
}

test <- test[!(is.na(test$from) & is.na(test$to)),]
# test <- test[complete.cases(test),]

g <- graph_from_edgelist(as.matrix(test[-c(1)]), directed = T) %>% 
  set_edge_attr("weight", value = test$weight)

plot(g)
require(visNetwork)
data <- toVisNetworkData(g)
visNetwork(nodes = data$nodes, edges = data$edges) %>% 
  visEdges(arrows = "to")



nodes <- data.frame(id = c(na.omit(unique(test$from, test$to)), NA))
edges <- data.frame(from = test$from, to = test$to)
visNetwork(nodes, edges, width = "100%")


##############hmm
dat <- tibble::tribble(
  ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`, ~`11`, ~`12`,
  10L,  "A",  "A",  "A",  "G",   NA,   NA,   NA,   NA,   NA,    NA,    NA,    NA,
  90L,  "B",  "A",  "A",  "C",   NA,   NA,   NA,   NA,   NA,    NA,    NA,    NA,
  17L,   NA,   NA,   NA,  "G",  "A",  "A",   NA,   NA,   NA,    NA,    NA,    NA,
  83L,   NA,   NA,   NA,  "C",  "T",  "A",   NA,   NA,   NA,    NA,    NA,    NA,
  30L,   NA,   NA,   NA,   NA,   NA,   NA,  "A",  "C",  "–",   "G",    NA,    NA,
  70L,   NA,   NA,   NA,   NA,   NA,   NA,  "A",  "C",  "A",   "G",    NA,    NA,
  50L,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,   "G",   "G",
  50L,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,   "G",   "D"
)



require(HMM)
# Initialise HMM
hmm = initHMM(c("A","B"), c("L","R"), transProbs=matrix(c(.6,.4,.4,.6),2),
              emissionProbs=matrix(c(.6,.4,.4,.6),2))
print(hmm)
# Sequence of observations
observations = c("L","L","R","R")
# Calculate Viterbi path
viterbi = viterbi(hmm, observations)
print(viterbi)




# Initial HMM
hmm = initHMM(c("A","B"),c("L","R"),
              transProbs=matrix(c(.9,.1,.1,.9),2),
              emissionProbs=matrix(c(.5,.51,.5,.49),2))
print(hmm)
# Sequence of observation
a = sample(c(rep("L",100),rep("R",300)))
b = sample(c(rep("L",300),rep("R",100)))
observation = c(a,b)
# Viterbi-training
vt = viterbiTraining(hmm,observation,10)
print(vt$hmm)






simHMM





tmp <- dat[-c(1)]
tmp.freq <- freq[-c(1)]


colnames(tmp) <- colnames(tmp.freq) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M")




hmm <- initHMM(States = colnames(tmp),
               Symbols = tmp,
               transProbs=matrix(tmp.freq))

tmp.obs <- ?
  (viterbi = viterbi(hmm, tmp.obs))
(vt = viterbiTraining(hmm,tmp.obs,10))



install.packages('DEploid')
























## Read example FASTA file.
require(haplotypes)
f<-system.file("example.fas",package="haplotypes")
# invalid character 'N' was replaced with '?' with a warning message
x<-read.fas(file=f)
# an object of class 'Dna'
x
## or load DNA Sequence data set.
data("dna.obj")
x<-dna.obj
## Not run:
x
## End(Not run)
## Compute an absolute pairwise character difference matrix from DNA sequences.
# coding gaps using simple indel coding method
d<- distance(x,indels="sic")
## Not run:
d
## End(Not run)
## Infer haplotypes using the 'Dna' object.
# coding gaps using simple indel coding method
h<-haplotype(x,indels="s")
## Not run:
h
## End(Not run)
## Conduct statistical parsimony analysis with 95% connection limit.
#algortihmic method
## Not run:
p<-parsimnet(x,prob=.95)
plot(p)








#############uninnnn

require(tidyverse)
dat <- tibble::tribble(
  ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`, ~`11`, ~`12`,
  10L,  "A",  "A",  "A",  "G",   NA,   NA,   NA,   NA,   NA,    NA,    NA,    NA,
  90L,  "B",  "A",  "A",  "B",   NA,   NA,   NA,   NA,   NA,    NA,    NA,    NA,
  17L,   NA,   NA,   NA,  "G",  "A",  "A",   NA,   NA,   NA,    NA,    NA,    NA,
  83L,   NA,   NA,   NA,  "B",  "T",  "A",   NA,   NA,   NA,    NA,    NA,    NA,
  30L,   NA,   NA,   NA,   NA,   NA,   NA,  "A",  "C",  "A",   "G",    NA,    NA,
  70L,   NA,   NA,   NA,   NA,   NA,   NA,  "A",  "C",  "A",   "B",    NA,    NA,
  50L,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,   "G",   "G",
  50L,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,   "G",   "D"
)

dat$`0` <- ifelse(is.na(dat$`1`), NA, "Y")
dat <- dat %>% 
  select(Count, `0`, everything())

freq <- data.frame(dat[c(1)])
for (c in 2:ncol(dat)){
  d <- dat[,c(1, c)]
  d.sum <- sum(d[complete.cases(d),]$Count)
  d$Count[is.na(d[c(2)])] <- 0
  d <- d$Count/d.sum
  freq <- cbind(freq, d)
}
colnames(freq) <- colnames(dat)


# if NA before or after then skip, elase, add _row
for (r in 1:nrow(dat)){
  for (c in 2:ncol(dat)){
    if (c < ncol(dat) -1){
      if (is.na(dat[r, c-1]) | is.na(dat[r, c+1])){
        #arg
      }else{
        dat[r, c] <- paste0(dat[r, c], "_", r)
      }
    }else{
      dat[r, c] <- paste0(dat[r, c], "_", r)
    }
  }
}

for (c in 2:ncol(dat)){
  dat[,c] <- paste0(unlist(dat[,c]), "_", colnames(dat)[c])
}
for (r in 1:nrow(dat)){
  dat[r,][grepl("NA_", dat[r,])] <- NA
}


test <- NULL
for (c in seq(ncol(dat)-1)){
  if (c>1){
    test.tmp <- dat[,c(1, c, c+1)]
    test.tmp$Count <- unlist(freq[c])
    colnames(test.tmp) <- c("weight", "from", "to")
    test <- rbind(test, test.tmp)
  }
}


test <- test[!(is.na(test$from)),] #test <- test[complete.cases(test),]

# if there is a from to a NA, and that from also occurs to a non NA to, then add an equal portion to those occurences
# i dont think thsi shoud be an equal proportion, come fix
tmp.join <- test[is.na(test$to),]
for (i in unique(tmp.join$from)){
  # i=unique(tmp.join$from)[1]
  if (sum(i==test$to, na.rm = T) == 1){
    # print(i)
    # print(sum(test[test$from == i, "weight"]))
    test[test$from == i, "weight"] <- sum(test[test$from == i, "weight"])
  }
}
tmp.join <- tmp.join[tmp.join$from%in%test$to,]
test <- test[complete.cases(test),]
# add  first, Y, and final, Z connection
all_nodes <- unique(c(test$from, test$to))
missing_conn_out <- data.frame()
# add in missing connections
missing_conn <- all_nodes[!all_nodes%in%test$from]
missing_conn <- missing_conn[!grepl("_0$", missing_conn)]
# missing_conn <- missing_conn[-length(missing_conn)]


for (node in missing_conn){
  pos <- as.integer(tail(str_split(node, "_")[[1]], 1))
  if (pos == as.integer(colnames(dat)[ncol(dat)])){
    nxt_pos <- "Z"
  }else{
    nxt_pos <- all_nodes[grepl(paste0("_", pos+1, "$"), all_nodes)]
  }
  for (to in nxt_pos){
    if (length(nxt_pos)==1){
      missing_conn_out <- rbind(missing_conn_out, 
                                c(1, node, to))
    }else{
      # get the sum of the connections from the next node
      to_weight <- sum(test[grepl(paste0(to, "$"), test$from), "weight"])
      missing_conn_out <- rbind(missing_conn_out, 
                                c(to_weight, node, to))
    }
  }
}

colnames(missing_conn_out) <- colnames(test)
missing_conn_out$weight <- as.numeric(missing_conn_out$weight)
test <- rbind(test,
              missing_conn_out)


test <- test %>% 
  group_by(from) %>% 
  mutate(weight = ifelse(weight>=1, weight/sum(weight), weight))




test$from[grepl("^Y", test$from)] <- "Y_0"
require(igraph)
g <- graph_from_edgelist(as.matrix(test[-c(1)]), directed = T) %>% 
  set_edge_attr("weight", value = test$weight)

# plot(g)
require(visNetwork)
data <- toVisNetworkData(g)
data$edges$label <- round(data$edges$weight, 3)
data$edges$value <- data$edges$weight*5
data$nodes$color.highlight.border <- "red"
n <- visNetwork(nodes = data$nodes, edges = data$edges) %>%
  visEdges(arrows = "to") %>% 
  visPhysics(enabled = F)
n %>% 
  visInteraction(multiselect=T, selectable=T) %>% 
  visOptions(selectedBy = list(
    variable = "label",
    style = "width:500px",
    multiple = TRUE, sort = FALSE
  ))

run <- 0
out_paths <- out_min_weight <- list()
out_remainder <- NULL
delete_this <- NULL
delete_this_too <- g
while (run <= length(all_simple_paths(g, from = "Y_0", to = "Z"))){
  g.tmp <- g
  # make largest smallest
  E(g.tmp)$weight <- 1/(E(g.tmp)$weight) #/max(E(g)$weight))
  p <- shortest_paths(g.tmp, from = "Y_0", to = "Z", mode = "out")$vpath[[1]]
  p.weights <- E(g, path = p)$weight
  # n.nodes <- length(p.weights)
  p.weights <- p.weights #[-n.nodes]
  if (!any(p.weights<=0)){
    out_paths <- out_paths %>% 
      append(list(p))
    out_min_weight <- out_min_weight %>% 
      append(min(p.weights))
    
    
    p.weights.delete <- E(delete_this_too, path = p)$weight
    delete_this <- c(delete_this,
                     prod((p.weights.delete/14))
                     )
    
    
    
    # update weights
    E(g, path = p)$weight <- p.weights-min(p.weights)
    E(g, path = p)$weight <- ifelse(E(g, path = p)$weight < 0, 0, E(g, path = p)$weight)
    # print(E(g, path = p)$weight)
  }else{
    run <- run + 1
  }
}

out_remainder <- g

for (i in 1:length(out_paths)){
  out_paths[[i]] <- names(out_paths[[i]])
}
df_out <- cbind(freq = unlist(out_min_weight), 
                delete_this = delete_this,
           do.call(rbind.data.frame, out_paths))

# print(E(g)$weight)
# out_paths
# out_min_weight

data2 <- toVisNetworkData(g)
data2$edges$label <- round(data2$edges$weight, 3)
data2$edges$value <- data2$edges$weight*5
data2$nodes$color.highlight.border <- "red"
n2 <- visNetwork(nodes = data2$nodes, edges = data2$edges) %>%
  visEdges(arrows = "to") %>% 
  visPhysics(enabled = F) %>% 
  visInteraction(multiselect=T, selectable=T)

n2












#############################uniteled next










haps_out <- data.frame()

for (chr in unique(phase_tmp$Chrom)){
  # chr="rpoB"
  
  
  split_in <- dat$Seq %>%
    str_split_fixed("", nchar(.))
  
  dat <- cbind(dat[-c(1)], split_in)
  dat[dat == "N"] <- NA
  
  # if NA before or after then skip, elase, add _row
  for (r in 1:nrow(dat)){
    for (c in 2:ncol(dat)){
      if (c < ncol(dat) -1){
        if (is.na(dat[r, c-1]) | is.na(dat[r, c+1])){
          #arg
        }else{
          dat[r, c] <- paste0(dat[r, c], "_", r)
        }
      }else{
        dat[r, c] <- paste0(dat[r, c], "_", r)
      }
    }
  }
  
  for (c in 2:ncol(dat)){
    dat[,c] <- paste0(unlist(dat[,c]), "_", colnames(dat)[c])
  }
  for (r in 1:nrow(dat)){
    dat[r,][grepl("NA_", dat[r,])] <- NA
  }
  
  test <- NULL
  for (c in seq(ncol(dat)-1)){
    if (c>1){
      test.tmp <- dat[,c(1, c, c+1)]
      colnames(test.tmp) <- c("weight", "from", "to")
      test <- rbind(test, test.tmp)
    }
  }
  
  test <- test[!(is.na(test$from) & is.na(test$to)),]
  test <- test[complete.cases(test),]
  
  test <- test %>% 
    group_by(from, to) %>% 
    mutate(weight = sum(weight)) %>% 
    unique()
  
  
  
  
  # add  first, Y, and final, Z connection
  all_nodes <- unique(c(test$from, test$to))
  missing_conn_out <- data.frame()
  
  for (first in all_nodes[grepl("_1$", all_nodes)]){
    missing_conn_out <- rbind(missing_conn_out, 
                              c(1, "Y", first))
  }
  
  # add in missing connections
  missing_conn <- all_nodes[!all_nodes%in%test$from]
  missing_conn <- missing_conn[!grepl("_1$", missing_conn)]
  # missing_conn <- missing_conn[-length(missing_conn)]
  
  
  for (node in missing_conn){
    pos <- as.integer(tail(str_split(node, "_")[[1]], 1))
    if (pos == as.integer(colnames(dat)[ncol(dat)])){
      nxt_pos <- "Z"
    }else{
      nxt_pos <- all_nodes[grepl(paste0("_", pos+1, "$"), all_nodes)]
    }
    for (to in nxt_pos){
      missing_conn_out <- rbind(missing_conn_out, 
                                c(1, node, to))
    }
  }
  colnames(missing_conn_out) <- colnames(test)
  missing_conn_out$weight <- as.integer(missing_conn_out$weight)
  test <- rbind(test,
                missing_conn_out)
  
  
  # test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1])
  
  g <- graph_from_edgelist(as.matrix(test[-c(1)]), directed = T) %>% 
    set_edge_attr("weight", value = test$weight)
  
  # plot(g)
  require(visNetwork)
  data <- toVisNetworkData(g)
  data$edges$label <- data$edges$weight
  # visNetwork(nodes = data$nodes, edges = data$edges) %>% 
  #   visEdges(arrows = "to")
  
  
  
  
  ps <- all_simple_paths(g, from = "Y", to = "Z")
  dat2 <- NULL
  # remove test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1]) is using this
  for (p in ps){
    dat2 <- rbind(dat2,
                  c(Count = sum(strength(g, vids = p, mode = "out")[-1]),
                    gsub("_.*", "", names(p))[-1])
    )
  }
  dat2 <- data.frame(dat2)
  
  # 
  # 
  # 
  # 
  # 
  # 
  # dat <- dat2
  # dat$Count <- as.integer(dat$Count)
  # r <- dat2[-c(1)]
  # frequency_assignment <- matrix(0, nrow = nrow(r), ncol = nrow(dat))
  # for(x in 1:nrow(dat)){
  #   valid <- apply(r, 1, function(z){
  #     all(z == dat[x, -1], na.rm = T)
  #   })
  #   frequency_assignment[valid , x] = dat$Count[x]/sum(valid)
  # }
  # frequency_assignment
  # 
  # sequences <- data.frame(frequency = rowSums(frequency_assignment),
  #                         probability = rowSums(frequency_assignment)/sum(frequency_assignment),
  #                         r)
  # sequences
  # 
  
  
  
  # walk through graph, taking next step that minimizus the edge - next edge
  out <- NULL
  for (p in ps){
    s <- strength(g, vids = p, mode = "out")[-1]
    ot <- 0
    for (i in 1:length(s)){
      v <- ifelse(s[i] == 1, s[i-1], s[i])
      ot <- v - ot
    }
    out <- c(out, ot)
  }
  # out
  
  
  dat2$cost <- out
  dat2$Count <- as.integer(dat2$Count)
  dat2$Count <- dat2$Count/sum(dat2$Count)
  # dat2
  
  # dat2 <- dat2 %>% 
  #   # select_if(~n_distinct(.) > 1) %>% 
  #   arrange(cost)
  # dat2
  
  
  # we want this many graphs
  need <- max(degree(g)-1)
  # issue with this is duplicates.... so fixed in the second dat 2 below
  # and had to add the if last column (c == ncol(dat2)-1)
  
  
  
  
  
  
  
  
  
  
  for (c in 2:ncol(dat2)-1){
    if (n_distinct(dat2[,c]) == need | (c == ncol(dat2)-1)){
      dat2 <- dat2 %>% 
        group_by_at(1:c) %>% 
        filter(cost == min(cost, na.rm = T)) %>% #whay are there NAs here?
        # slice(1) %>% # takes the first occurrence if there is a tie
        ungroup() %>% 
        # mutate(Count = Count/sum(Count)) %>% 
        select(cost, everything())
      break
    }
  }
  
  
  dat2 <- dat2 %>% 
    group_by_at(1:c) %>% 
    mutate(split = dplyr::cur_group_id()) %>%   #need to mark ties, so they are not counted as unique for the Count
    filter(cost == min(cost, na.rm = T)) %>% 
    unique() # this added as code before is bad, very dangerous!
  
  seqs <- dat2 %>% 
    ungroup() %>% 
    select(-cost, -Count, -split, -ncol(.)) %>% 
    select(-ncol(.)) %>% 
    unite("Seq", remove = F, sep = "") %>% 
    select(Seq)
  
  dat2 <- dat2 %>% 
    ungroup() %>% 
    select(cost, Count, split) %>% 
    mutate(Seq = unlist(seqs))
  
  dat2$split <- ifelse(duplicated(dat2$split), dat2$split, "one")
  
  dat2 <- dat2 %>% 
    group_by(Seq, split) %>% 
    mutate(Count = sum(Count)) %>% 
    filter(cost == min(cost, na.rm = T)) %>% 
    ungroup()
  
  dat2 <- dat2 %>% 
    mutate(Count = Count/sum(Count)) %>% 
    unique()
  dat2 <- dat2 %>% 
    mutate(Count = Count/sum(Count)) %>% 
    unique()
  
  
  
  haps_out <- rbind(haps_out, 
                    cbind(chr, dat2))
  
}








#####################untieled nnn

dat <- tibble::tribble(
  ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`,~`11`,~`12`,
  1001L,  "A",  "A",  "A",  "G",  NA,  NA,  NA,  NA,  NA,   NA,  NA,NA,
  506L,  NA,  NA,  NA,  "G",  "A",  "A",  NA,  NA,  NA,   NA, NA,NA,
  1030L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "–",   "G", NA,NA,
  1030L,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, "G", "G",
  
  9020L,  "B",  "A",  "A",  "C",  NA,  NA,  NA,  NA,  NA,   NA,  NA,NA,
  4035L,  NA,  NA,  NA,  "C",  "T",  "A",  NA,  NA,  NA,   NA, NA,NA,
  4025L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "A",   "G", NA,NA,
  4025L,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, "G", "D"
)




# dat <- dat[apply(dat, 2, function(x) length(unique(na.omit(x))))>1]
# colnames(dat) <- c("Count", 1:ncol(dat))

# make ratios by read length
# dat$grp <- "_"
# for (r in 1:nrow(dat)){
#   first_na <- min(which(is.na(dat[r,])))
#   last_na <- max(which(is.na(dat[r,])))
#   dat[r,"grp"] <- paste0(first_na, "_", last_na)
# }

# dat <- dat %>%
#   group_by(grp) %>%
#   mutate(Count = Count/sum(Count)) %>%
#   ungroup()
# 
# dat <- dat %>%
#   select(-grp)

dat$`0` <- ifelse(is.na(dat$`1`), NA, "Y")
dat <- dat %>% 
  select(Count, `0`, everything())

freq <- data.frame(dat[c(1)])
for (c in 2:ncol(dat)){
  d <- dat[,c(1, c)]
  d.sum <- sum(d[complete.cases(d),]$Count)
  d$Count[is.na(d[c(2)])] <- 0
  d <- d$Count/d.sum
  freq <- cbind(freq, d)
}
colnames(freq) <- colnames(dat)







# if NA before or after then skip, elase, add _row
for (r in 1:nrow(dat)){
  for (c in 2:ncol(dat)){
    if (c < ncol(dat) -1){
      if (is.na(dat[r, c-1]) | is.na(dat[r, c+1])){
        #arg
      }else{
        dat[r, c] <- paste0(dat[r, c], "_", r)
      }
    }else{
      dat[r, c] <- paste0(dat[r, c], "_", r)
    }
  }
}

for (c in 2:ncol(dat)){
  dat[,c] <- paste0(unlist(dat[,c]), "_", colnames(dat)[c])
}
for (r in 1:nrow(dat)){
  dat[r,][grepl("NA_", dat[r,])] <- NA
}


test <- NULL
for (c in seq(ncol(dat)-1)){
  if (c>1){
    test.tmp <- dat[,c(1, c, c+1)]
    test.tmp$Count <- unlist(freq[c])
    colnames(test.tmp) <- c("weight", "from", "to")
    test <- rbind(test, test.tmp)
  }
}


test <- test[!(is.na(test$from)),] #test <- test[complete.cases(test),]

# 
# test <- test %>%
#   group_by(from) %>%
#   mutate(weight = sum(weight)/n()) %>%
#   unique()
# 
# 
# test$PosFrom <- 0
# for (r in 1:nrow(test)){
#   test[r, "PosFrom"] <- as.integer(tail(str_split(test[r,]$from, "_")[[1]], 1))
# }
# 
# test <- test %>%
#   group_by(from) %>%
#   mutate(weight = sum(weight))

test <- test[complete.cases(test),]





# add  first, Y, and final, Z connection
all_nodes <- unique(c(test$from, test$to))
missing_conn_out <- data.frame()

# for (first in all_nodes[grepl("_1$", all_nodes)]){
#   missing_conn_out <- rbind(missing_conn_out, 
#                             c(1, "Y", first))
# }

# add in missing connections
missing_conn <- all_nodes[!all_nodes%in%test$from]
missing_conn <- missing_conn[!grepl("_0$", missing_conn)]
# missing_conn <- missing_conn[-length(missing_conn)]


for (node in missing_conn){
  pos <- as.integer(tail(str_split(node, "_")[[1]], 1))
  if (pos == as.integer(colnames(dat)[ncol(dat)])){
    nxt_pos <- "Z"
  }else{
    nxt_pos <- all_nodes[grepl(paste0("_", pos+1, "$"), all_nodes)]
  }
  for (to in nxt_pos){
    
    val_prvious <- 
      
      missing_conn_out <- rbind(missing_conn_out, 
                                c(1, node, to))
  }
}
colnames(missing_conn_out) <- colnames(test)
missing_conn_out$weight <- as.integer(missing_conn_out$weight)
test <- rbind(test,
              missing_conn_out)






test$from[grepl("^Y", test$from)] <- "Y_0"
# if only one in then make it 1
test <- test %>%
  group_by(to) %>%
  mutate(weight2 = ifelse(n() == 1 , 1, weight)) %>% #& from != "Y_0"
  group_by(to, from) %>% 
  mutate(weight = ifelse(weight2 == 1, 1, weight)) %>% 
  ungroup() %>% 
  select(-weight2)







# test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1])

g <- graph_from_edgelist(as.matrix(test[-c(1)]), directed = T) %>% 
  set_edge_attr("weight", value = test$weight)

# plot(g)
require(visNetwork)
data <- toVisNetworkData(g)
data$edges$label <- round(data$edges$weight, 3)
visNetwork(nodes = data$nodes, edges = data$edges) %>%
  visEdges(arrows = "to") %>% 
  visPhysics(enabled = F)





ps <- all_simple_paths(g, from = "Y_0", to = "Z")
dat2 <- NULL
# remove test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1]) is using this
for (p in ps){
  dat2 <- rbind(dat2,
                c(Count = sum(strength(g, vids = p, mode = "out")[-1]),
                  gsub("_.*", "", names(p))[-1])
  )
}
dat2 <- data.frame(dat2)

dat2$Count <- as.numeric(dat2$Count)/(length(ps[[1]])-3)


# dat2$Count <- dat2$Count/sum(dat2$Count)

# 
# dat2
# freq
# 
# dim(dat2[-c(ncol(dat2))])
# dim(freq)
# 


out <- NULL
for (p in ps){
  # s <- sort(strength(g, vids = p, mode = "out"))[2]
  s1 <- strength(g, vids = p, mode = "out")[-c(1, length(p))]
  s1[s1>1] <- 1 # fix this in the code above
  s <- prod(sort(
    s1
  )) #/length(p)
  out <- c(out, s)
}


dat2$cost <- out
# dat2$Count <- as.integer(dat2$Count)
# dat2$Count <- dat2$Count/sum(dat2$Count)
dat2





# 
# # walk through graph, taking next step that minimizus the edge - next edge
# out <- NULL
# for (p in ps){
#   s <- strength(g, vids = p, mode = "out")[-1]
#   ot <- 0
#   for (i in 2:length(s)){
#     if (s[i] != 1 & s[i-1] != 1){
#       ot <- ot + (s[i-1]%%s[i])
#     }
#   }
#   out <- c(out, ot)
# }
# out

# 
# 
# 
# s <- strength(g, vids = p, mode = "out")[-1]
# cst <- 1
# for (i in s){
#   if (i != 0 & i != 1){
#     if (cst == 0){
#       cst <- i
#     }else{
#       print(paste("i", i))
#       print(paste("rem", cst %% i))
#       cst <- cst + (cst %% i)
#       print(paste("cst", cst))
#     }
#   }
#   cst <- cst - 1
# }
# s




##################uniteles nnn

dat <- tibble::tribble(
  ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`,~`11`,~`12`,
  1001L,  "A",  "A",  "A",  "G",  NA,  NA,  NA,  NA,  NA,   NA,  NA,NA,
  506L,  NA,  NA,  NA,  "G",  "A",  "A",  NA,  NA,  NA,   NA, NA,NA,
  1030L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "–",   "G", NA,NA,
  1030L,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, "G", "G",
  
  9020L,  "B",  "A",  "A",  "C",  NA,  NA,  NA,  NA,  NA,   NA,  NA,NA,
  4035L,  NA,  NA,  NA,  "C",  "T",  "A",  NA,  NA,  NA,   NA, NA,NA,
  4025L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "A",   "G", NA,NA,
  4025L,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, "G", "D"
)




# dat <- dat[apply(dat, 2, function(x) length(unique(na.omit(x))))>1]
# colnames(dat) <- c("Count", 1:ncol(dat))

# make ratios by read length
# dat$grp <- "_"
# for (r in 1:nrow(dat)){
#     first_na <- min(which(is.na(dat[r,])))
#     last_na <- max(which(is.na(dat[r,])))
#     dat[r,"grp"] <- paste0(first_na, "_", last_na)
# }
# 
# dat <- dat %>% 
#   group_by(grp) %>% 
#   mutate(Count = Count/sum(Count)) %>% 
#   ungroup()

# dat <- dat %>% 
#   select(-grp)
freq <- data.frame(dat[c(1)])
for (c in 2:ncol(dat)){
  d <- dat[,c(1, c)]
  d.sum <- sum(d[complete.cases(d),]$Count)
  d$Count[is.na(d[c(2)])] <- 0
  d <- d$Count/d.sum
  freq <- cbind(freq, d)
}
colnames(freq) <- colnames(dat)



# if NA before or after then skip, elase, add _row
for (r in 1:nrow(dat)){
  for (c in 2:ncol(dat)){
    if (c < ncol(dat) -1){
      if (is.na(dat[r, c-1]) | is.na(dat[r, c+1])){
        #arg
      }else{
        dat[r, c] <- paste0(dat[r, c], "_", r)
      }
    }else{
      dat[r, c] <- paste0(dat[r, c], "_", r)
    }
  }
}

for (c in 2:ncol(dat)){
  dat[,c] <- paste0(unlist(dat[,c]), "_", colnames(dat)[c])
}
for (r in 1:nrow(dat)){
  dat[r,][grepl("NA_", dat[r,])] <- NA
}


test <- NULL
for (c in seq(ncol(dat)-1)){
  if (c>1){
    test.tmp <- dat[,c(1, c, c+1)]
    test.tmp$Count <- unlist(freq[c])
    colnames(test.tmp) <- c("weight", "from", "to")
    test <- rbind(test, test.tmp)
  }
}

test <- test[!(is.na(test$from) & is.na(test$to)),]
test <- test[complete.cases(test),]

test <- test %>% 
  group_by(from, to) %>% 
  mutate(weight = sum(weight)) %>% 
  unique()




# add  first, Y, and final, Z connection
all_nodes <- unique(c(test$from, test$to))
missing_conn_out <- data.frame()

for (first in all_nodes[grepl("_1$", all_nodes)]){
  missing_conn_out <- rbind(missing_conn_out, 
                            c(1, "Y", first))
}

# add in missing connections
missing_conn <- all_nodes[!all_nodes%in%test$from]
missing_conn <- missing_conn[!grepl("_1$", missing_conn)]
# missing_conn <- missing_conn[-length(missing_conn)]


for (node in missing_conn){
  pos <- as.integer(tail(str_split(node, "_")[[1]], 1))
  if (pos == as.integer(colnames(dat)[ncol(dat)])){
    nxt_pos <- "Z"
  }else{
    nxt_pos <- all_nodes[grepl(paste0("_", pos+1, "$"), all_nodes)]
  }
  for (to in nxt_pos){
    missing_conn_out <- rbind(missing_conn_out, 
                              c(1, node, to))
  }
}
colnames(missing_conn_out) <- colnames(test)
missing_conn_out$weight <- as.integer(missing_conn_out$weight)
test <- rbind(test,
              missing_conn_out)


# test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1])

g <- graph_from_edgelist(as.matrix(test[-c(1)]), directed = T) %>% 
  set_edge_attr("weight", value = test$weight)

# plot(g)
require(visNetwork)
data <- toVisNetworkData(g)
data$edges$label <- round(data$edges$weight, 2)
visNetwork(nodes = data$nodes, edges = data$edges) %>%
  visEdges(arrows = "to") %>% 
  visPhysics(enabled = F)





ps <- all_simple_paths(g, from = "Y", to = "Z")
dat2 <- NULL
# remove test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1]) is using this
for (p in ps){
  dat2 <- rbind(dat2,
                c(Count = sum(strength(g, vids = p, mode = "out")[-1]),
                  gsub("_.*", "", names(p))[-1])
  )
}
dat2 <- data.frame(dat2)

dat2$Count <- as.numeric(dat2$Count)/(length(ps[[1]])-3)


# dat2$Count <- dat2$Count/sum(dat2$Count)


dat2
freq

dim(dat2[-c(ncol(dat2))])
dim(freq)
