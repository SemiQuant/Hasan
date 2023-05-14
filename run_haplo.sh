echo bam='"'${1}'"' > inputs.R
echo snps='"'${2}'"' >> inputs.R

Rscript -e "rmarkdown::render('phasing_haplo.Rmd')"