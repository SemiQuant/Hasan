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

form this i want the most likley full sequences (rows) and frequencies, based on the counts and the sequences, here it would be
dat_out <- tibble::tribble(
             ~Frequency, ~`1`, ~`2`, ~`3`, ~`4`, ~`5`, ~`6`, ~`7`, ~`8`, ~`9`, ~`10`,
                    0.1,  "A",  "A",  "A",  "G",  "A",  "A",  "A",  "C",  "–",   "G",
                    0.9,  "A",  "A",  "A",  "C",  "T",  "A",  "A",  "C",  "A",   "G"
             )

It wont always be that easy, but i think this is a pretty generalizable example
(like the non overlapping seuqneces may be at the same frequencey and thus both options will be 'just as likley')
And sometime there will be 
AA'NA'A
and 
AAAA
etc.

something like this
https://ars.els-cdn.com/content/image/1-s2.0-S0168170216304130-gr3.jpg