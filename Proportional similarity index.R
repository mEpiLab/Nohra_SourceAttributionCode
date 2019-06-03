library(dplyr)
library(stringr)
#library("devtools")
#install_github("hadley/dplyr")

# read data in
#setwd("H:/Antoine- Massey/Thesis/C. jejuni/PSI_C. jejuni")
setwd ("C:/Users/Antoine/Documents/Antoine- Massey/Thesis/C. jejuni/PSI_C. jejuni")

db <- read.csv("all_sources_after_manawatu.csv", stringsAsFactors=F)

# change this to the source mapping file you want
source_file <- "all_source"

#sources <- read.csv(file.path("H:/Antoine- Massey/Thesis/C. jejuni/PSI_C. jejuni/input", paste0(source_file, ".csv")), colClasses="character")

sources <- read.csv(file.path("C:/Users/Antoine/Documents/Antoine- Massey/Thesis/C. jejuni/PSI_C. jejuni/input", paste0(source_file, ".csv")), colClasses="character")
source_map <- as.numeric(sources$Number)
names(source_map) <- sources$DataSource

source_label_map <- unique(sources %>% select(Number, Label))
source_labels <- str_replace(source_label_map$Label, "\\\\n", "\n")
names(source_labels) <- source_label_map$Number

humans <- db %>% filter(Source == "Human")
animals <- db %>% filter(Source %in% sources$DataSource)
animals <- animals %>% mutate(Source = source_labels[as.character(source_map[as.character(Source)])])

db <- rbind(humans, animals)

tab <- table(db$Source, db$ST)

# table of ST counts per source
options(max.print = 5.5E5)
print(tab)

# convert table to frequencies (proportions)
ptab <- prop.table(tab, margin=1)

# work out the PSI
psi <- matrix(NA, nrow(ptab), nrow(ptab))

for (i in 1:nrow(ptab)) {
  for (j in 1:nrow(ptab)) {
    psi[i,j] = 1 - 0.5*sum(abs(ptab[i,] - ptab[j,]))
  }
}

rownames(psi) <- rownames(ptab)
colnames(psi) <- rownames(ptab)


####################Similarity index Confidence intervals############################
#setwd("H:/Antoine- Massey/Thesis/C. jejuni/PSI_C. jejuni")

setwd ("C:/Users/Antoine/Documents/Antoine- Massey/Thesis/C. jejuni/PSI_C. jejuni")

similarityIndex<-function(x,y) {1-0.5*sum(abs(x/sum(x)-y/sum(y)))}
getSample<-function(x) {
  sa<-sample(c(1:length(x)),sum(x),replace=TRUE,prob=x/sum(x))
  out<-matrix(0,1,length(x))
  for (i in 1:length(x)) {
    out[i]<-length(sa[sa==i])
  }
  out
}

iters<-2000
sources<-7   # Put the number of sources you have
# COlumns should be the sources without adding them in the excel file."Humans, Poultry"
# Rows are the ST types and also should not be added. Just add the frequency of each ST type for each source


# Load data - ST assigned for each source

#d<-scan('all_sources_after.txt')  
d<-scan('all_sources_after_manawatu.txt')

types<-length(d)/sources
isolates<-matrix(d,sources,types)
combinations<-sources*(sources-1)/2
comb<-matrix(0,combinations,2)
k<-1
for (i in 1:(sources-1)) {
  for (j in (i+1):sources) {
    comb[k,]<-c(i,j)
    k<-k+1
  }
}
simin<-matrix(0,combinations,iters)
for (i in 1:combinations) {
  for (j in 1:iters) {
    simin[i,j]<-similarityIndex(getSample(isolates[comb[i,1],]),getSample(isolates[comb[i,2],]))
  }
  cat(comb[i,1]," <-> ",comb[i,2]," PS=",similarityIndex(isolates[comb[i,1],],isolates[comb[i,2],])," (",quantile(simin[i,],0.025),",",quantile(simin[i,],0.975),")\n",sep="")
}
