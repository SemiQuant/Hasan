dat <- tibble::tribble(
  ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`,
  1001L,  "A",  "A",  "A",  "G",  NA,  NA,  NA,  NA,  NA,   NA,
  9020L,  "A",  "A",  "A",  "C",  NA,  NA,  NA,  NA,  NA,   NA,
  506L,  NA,  NA,  NA,  "G",  "A",  "A",  NA,  NA,  NA,   NA,
  4035L,  NA,  NA,  NA,  "C",  "T",  "A",  NA,  NA,  NA,   NA,
  630L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "–",   "G",
  4025L,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "A",   "G"
)

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

plot(g)
require(visNetwork)
data <- toVisNetworkData(g)
data$edges$label <- data$edges$weight
visNetwork(nodes = data$nodes, edges = data$edges) %>% 
  visEdges(arrows = "to")




ps <- all_simple_paths(g, from = "Y", to = "Z")
dat2 <- NULL
# remove test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1]) is using this
for (p in ps){
  dat2 <- rbind(dat2,
                c(Count = sum(strength(g, vids = p, mode = "out")[-1]),
                  gsub("_.*", "", names(p))[-1])
  )
}
(dat2 <- data.frame(dat2))

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
out


dat2$cost <- out
dat2$Count <- as.integer(dat2$Count)
dat2$Count <- dat2$Count/sum(dat2$Count)
dat2

# dat2 <- dat2 %>% 
#   # select_if(~n_distinct(.) > 1) %>% 
#   arrange(cost)
# dat2

# we want this many graphs
need <- max(degree(g)-1)

for (c in 2:ncol(dat2)){
  if (n_distinct(dat2[,c]) == need){
    dat2 <- dat2 %>% 
      group_by(dat2[,c]) %>% 
      filter(cost == min(cost)) %>% 
      # slice(1) %>% # takes the first occurrence if there is a tie
      ungroup() %>% 
      mutate(Count = Count/sum(Count)) %>% 
      select(-ncol(.)) %>% 
      select(cost, everything())
    break
  }
}

dat2

