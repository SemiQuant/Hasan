require(tidyverse)
require(igraph)
require(visNetwork)

# dat <- tibble::tribble(
#   ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`,~`11`,~`12`,
#   10,  "A",  "A",  "A",  "G",  NA,  NA,  NA,  NA,  NA,   NA,  NA,NA,
#   10,  NA,  NA,  NA,  "G",  "A",  "A",  NA,  NA,  NA,   NA, NA,NA,
#   10,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "–",   "G", NA,NA,
#   10,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, "G", "G",
#   
#   90,  "B",  "A",  "A",  "C",  NA,  NA,  NA,  NA,  NA,   NA,  NA,NA,
#   90,  NA,  NA,  NA,  "C",  "T",  "A",  NA,  NA,  NA,   NA, NA,NA,
#   90,  NA,  NA,  NA,  NA,  NA,  NA,  "A",  "C",  "A",   "G", NA,NA,
#   90,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, "G", "D"
# )


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


dat <- tibble::tribble(
  ~Count, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`, ~`11`, ~`12`,
  10L,  "A",  "A",  "A",  "G",   NA,   NA,   NA,   NA,   NA,    NA,    NA,    NA,
  17L,   NA,   NA,   NA,  "G",  "A",  "A",   NA,   NA,   NA,    NA,    NA,    NA,
  30L,   NA,   NA,   NA,   NA,   NA,   NA,  "A",  "C",  "–",   "G",    NA,    NA,
  
  90L,  "B",  "A",  "A",  "C",   NA,   NA,   NA,   NA,   NA,    NA,    NA,    NA,
  83L,   NA,   NA,   NA,  "C",  "T",  "A",   NA,   NA,   NA,    NA,    NA,    NA,
  70L,   NA,   NA,   NA,   NA,   NA,   NA,  "A",  "C",  "A",   "G",    NA,    NA,
  
  50L,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,   "G",   "G",
  50L,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,   "G",   "D"
)

dat[dat == "N"] <- NA
dat


























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
tmp.join


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
# # if only one in then make it 1
# test <- test %>%
#   group_by(to) %>%
#   mutate(weight2 = ifelse(n() == 1 , 1, weight)) %>% #& from != "Y_0"
#   group_by(to, from) %>% 
#   mutate(weight = ifelse(weight2 == 1, 1, weight)) %>% 
#   ungroup() %>% 
#   select(-weight2)







# test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1])

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



# 
# # n %>% 
#   # visSave(file = "/Users/SemiQuant/Downloads/net1.html")
# 
# 
# ps <- all_simple_paths(g, from = "Y_0", to = "Z")
# dat2 <- NULL
# # remove test$weight <- 1/(test$weight/sort(unique(test$weight))[length(unique(test$weight))-1]) is using this
# for (p in ps){
#   dat2 <- rbind(dat2,
#                 c(Count = sum(strength(g, vids = p, mode = "out")[-1]),
#                   gsub("_.*", "", names(p))[-1])
#   )
# }
# dat2 <- data.frame(dat2)
# 
# dat2$Count <- as.numeric(dat2$Count)/(length(ps[[1]])-3)
# 
# 
# # dat2$Count <- dat2$Count/sum(dat2$Count)
# 
# # 
# # dat2
# # freq
# # 
# # dim(dat2[-c(ncol(dat2))])
# # dim(freq)
# # 
# 
# 
# out <- NULL
# for (p in ps){
#   # s <- sort(strength(g, vids = p, mode = "out"))[2]
#   s1 <- strength(g, vids = p, mode = "out")[-c(1, length(p))]
#   s1[s1>1] <- 1 # fix this in the code above
#   s <- prod(sort(
#     s1
#   )) #/length(p)
#   out <- c(out, s)
# }
# 
# 
# dat2$cost <- out
# # dat2$Count <- as.integer(dat2$Count)
# # dat2$Count <- dat2$Count/sum(dat2$Count)
# dat2
# 
# 
# 
# 
# 
# # 
# # # walk through graph, taking next step that minimizus the edge - next edge
# # out <- NULL
# # for (p in ps){
# #   s <- strength(g, vids = p, mode = "out")[-1]
# #   ot <- 0
# #   for (i in 2:length(s)){
# #     if (s[i] != 1 & s[i-1] != 1){
# #       ot <- ot + (s[i-1]%%s[i])
# #     }
# #   }
# #   out <- c(out, ot)
# # }
# # out
# 
# # 
# # 
# # 
# # s <- strength(g, vids = p, mode = "out")[-1]
# # cst <- 1
# # for (i in s){
# #   if (i != 0 & i != 1){
# #     if (cst == 0){
# #       cst <- i
# #     }else{
# #       print(paste("i", i))
# #       print(paste("rem", cst %% i))
# #       cst <- cst + (cst %% i)
# #       print(paste("cst", cst))
# #     }
# #   }
# #   cst <- cst - 1
# # }
# # s
# 
# 
# 
# 




# 
# # make largest smallest
# E(g)$weight <- 1/(E(g)$weight/max(E(g)$weight))
# 
# # for{
# # get shortest path, uses dijkstra for directed
# p <- shortest_paths(g, from = "Y_0", to = "Z", mode = "out")$vpath[[1]]
# 
# # update weights
# p.weights <- E(g, path = p)$weight
# E(g, path = p)$weight <- p.weights + min(p.weights)
# E(g, path = p)$weight <- ifelse(E(g, path = p)$weight < 0, 0, E(g, path = p)$weight)
# # }




# \\\change this to  ?
# max_flow(g, source = "Y_0", target = "Z", )

# flow from source to target is an assignment of non-negative real numbers to the edges of the graph, 
# satisfying two properties: (1) for each edge the flow (ie. the assigned number) is not more than the 
# capacity of the edge (the capacity parameter or edge attribute), (2) for every vertex, except the source 
# and the target the incoming flow is the same as the outgoing flow. The value of the flow is the incoming 
# flow of the target vertex. The maximum flow is the flow of maximum value.



# I think this will act to reweight those tht cross over to the min
test.tmp <- test
test.tmp$capacity <- test.tmp$weight
g.tmp <- g
E(g.tmp)$capacity <- as.integer(E(g.tmp)$weight*1000)
mf <- max_flow(g.tmp, source = "Y_0", target = "Z")
E(g)$weight <- mf$flow/1000







# this is basically ba slightly edited, bad version of the max flow algorithm
run <- 0
out_paths <- out_min_weight <- list()
out_remainder <- NULL
delete_this <- NULL
delete_this_too <- g

# intermediate_gs <- list()
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
    
    # intermediate_gs <- intermediate_gs %>% 
    #   append(g)
    
  }else{
    run <- run + 1
  }
}

out_remainder <- g

for (i in 1:length(out_paths)){
  out_paths[[i]] <- names(out_paths[[i]])
}
df_out <- cbind(freq = unlist(out_min_weight), 
                # delete_this = delete_this,
                do.call(rbind.data.frame, out_paths))



df_out$seq <- apply(df_out[,3:(ncol(df_out)-1)], 2, function(x) gsub("_.*", "", x)) %>% 
  as.data.frame() %>% 
  unite("Seq", sep = "")


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

