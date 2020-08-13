# build a markov model
# given all the pssible sequences, what is the prob that they came from this model?




require(bnlearn)
dag <- as.bn(g)
arcs(dag)
# E(g)$weight # check this!
graphviz.plot(dag, layout = "fdp")
plot(g)
# 
test <- tibble::tribble(
            ~Y_0,  ~A_2_2,  ~A_2_2,  ~C_4,  ~C_4,  ~T_4_5,  ~`A_6`,  ~`A_7`, ~`C_6_8`,  ~`A_6_9`,  ~`G_10`,  ~`G_8_11`, ~`D_8_12`,
          "B_2_1", "A_2_2", "A_2_3", "C_4", "T_4_5", "A_6", "A_7", "C_6_8",  "A_6_9", "G_10", "G_8_11", "D_8_12", "Z",
          "A_1_1", "A_1_2", "A_1_3", "G_4", "A_3_5", "A_6", "A_7", "C_5_8",  "–_5_9", "G_10", "G_7_11", "G_7_12", "Z",
          "B_2_1", "A_2_2", "A_2_3", "C_4", "T_4_5", "A_6", "A_7", "C_6_8",  "A_6_9", "G_10", "G_8_11", "D_8_12", "Z"
          )
# 
# 
# bn.fit.barchart(bn.g$T, main = "Travel",
#                 xlab = "Pr(T | R,O)", ylab = "")


bn.fit(dag, data = test, method = "mle")










# discret case
dag <- model2network("[C][B][D|B:C]")
A.lv <- c("M", "F")



B.prob <- array(c(0.30, 0.70), 
                dim = 2,
                dimnames = list(D = A.lv))
C.prob <- array(c(0.30, 0.70), 
                dim = 2,
                dimnames = list(D = A.lv))
D.prob <- array(c(0.30, 0.70), 
                dim = 2,
                dimnames = list(D = A.lv))

cpt <- list(B = B.prob, C = C.prob, D = D.prob)
bn <- custom.fit(dag, cpt)

cpquery(bn, event = (S == "M") & (T == "car"), 
        evidence = (E == "high") & (R == "small"))







