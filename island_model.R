# code for running per-year island attribution

library(RColorBrewer)
library(stringr)
library(dplyr)
library("devtools")
install_github("hadley/dplyr")

setwd("H:/Antoine- Massey/Thesis/C. jejuni/attribution/attribution")

source("helpers.R")

# load in our dataset (relative to current working directory)
db <- read.csv("all_sources.csv")

# change this to the source mapping file you want
source_file <- "all_source"

sources <- read.csv(file.path("island/input", paste0(source_file, ".csv")), colClasses="character")
source_map <- as.numeric(sources$Number)
names(source_map) <- sources$DataSource

source_label_map <- unique(sources %>% select(Number, Label))
source_labels <- str_replace(source_label_map$Label, "\\\\n", "\n")
names(source_labels) <- source_label_map$Number

# input parameters
alleles_to_impute <- max(c(suppressWarnings(as.numeric(sources$Imputed)), 0), na.rm=T)

# model fitting control
seeds     <- c(5,7,11,13,17)
num_iters <- 20000
thinning  <- 50

human   <- "Human"
mlst_cols  <- c("ST", "ASP", "GLN", "GLT", "GLY", "PGM", "TKT", "UNC")

# Setup data
db <- db %>% filter(Imputed <= alleles_to_impute)

humans <- db %>% filter(Source == human)

animals <- db %>% filter(Source %in% sources$DataSource)
animals <- animals %>% mutate(Source = source_map[as.character(Source)])
animals <- animals %>% select(one_of(c(mlst_cols, "Source")))

datasets <- list()

datasets[[source_file]] <- humans

# make the temp directory
temp_dir <- file.path("H:/Antoine- Massey/Thesis/C. jejuni/attribution/attribution/island/temp", source_file)
dir.create(temp_dir, showWarnings=F, recursive=T)

for (year in 1:length(datasets)) {

  # make the year directory
  year_dir <- file.path(temp_dir, names(datasets)[year])
  dir.create(year_dir, showWarnings=F)

  # create human dataset (we use the same animal dataset each time)
  h_year <- datasets[[year]] %>% select(one_of(mlst_cols))
  h_year <- h_year %>% mutate(Source = 0)

  data <- rbind(h_year, animals)

  # save data file
  island_in <- "input.txt"
  write.table(data, file=file.path(year_dir, island_in), sep="\t", row.names=F)

  # copy our isource executable there (as isource is braindead about output shit)
  file.copy("island/isource_static.exe", year_dir)
  current_dir <- getwd()

  for (seed in seeds) {
    island_seed_out <- paste0("output_", seed, ".txt")
    command    <- paste0("isource_static ", island_in, " ", island_seed_out, " ", num_iters, " ", thinning, " 1 -", seed)

    cat("---------------------------------------\n")
    cat("running island model with seed", seed, "\n")
    cat(command, "\n")
    cat("---------------------------------------\n")

    setwd(year_dir)
    system(command)
    setwd(current_dir)
  }

  # analyse the output
  out_files <- list.files(year_dir, pattern = "^output_[0-9]+.txt")

  # TODO: improve so we don't have a dependence on animals$Source
  num_sources <- max(animals$Source)
  mcmc = NULL; fmcmc = NULL;
  for(i in seq_along(out_files)) {
    filename = out_files[i]
    mcmc = rbind(mcmc, read.table(file.path(year_dir, filename), header=T, comment.char=""))
    fmcmc_file = paste("f_", filename, sep="")
    fmcmc = rbind(fmcmc, read.table(file.path(year_dir, fmcmc_file), header=T, comment.char=""))
    g_file = paste("g_", filename, sep="")
    if (i==1) {
      g = matrix(scan(file.path(year_dir, g_file), what=double(0), sep="\t"), nrow = num_sources)
    } else {
      g = g + matrix(scan(file.path(year_dir, g_file), what=double(0), sep="\t"),nrow = num_sources)
    }
  }
  g = t(g)/length(out_files)

  # eliminate burnin
  gd = mcmc$iter>=5000
  fd = fmcmc$iter>=1000

  # produce table
  df_all = fmcmc[fd,2:(num_sources+1)]; names(df_all) <- NULL

  (pe = apply(df_all,2,function(x)c("mean"=mean(x),"median"=median(x),"sd"=sd(x),quantile(x,c(.025,.975)))))

  # write out the data
  data_out <- data.frame(source = 1:num_sources, mean = pe[1,], lci = pe[4,], uci = pe[5,], total=nrow(g))

  island_out <- file.path(year_dir, "output.txt")
  write.table(data_out, file=island_out, row.names=F)
}

# TODO: plotting, steal from island model single run...

# create output folder
output_dir <- file.path("H:/Antoine- Massey/Thesis/C. jejuni/attribution/attribution/island/output", "antoine")
dir.create(output_dir, showWarnings=F, recursive=T)

o <- c(1,5,4,6)
col = c("#FF7F00","#CF0000","#004FCF", "#009F9F","#8F006F","#9F5F3F","#FFAFAF")
col <- col[o]
plot_names <- names(datasets)

#plot_names <- c("Urban Rawmilk")

# function for creating an attribution plot
attribution_plot <- function(output_pdf, data_out, col=NULL, plot_name=NULL) {
  pdf(output_pdf, width=5, height=3.5)
  if (is.null(col)) {
    col <- rainbow(nrow(data_out))
  }
  mp = barplot(data_out$mean*100,col=col,ylim=c(0,100),ylab="Percentage of cases attributed", cex.axis=0.8, cex.lab=0.8)
  segments(mp,data_out$lci*100,mp,data_out$uci*100,lwd=2)

  for (i in 1:nrow(data_out)) {
    label <- source_labels[as.character(i)]
    carriage_returns <- nrow(str_locate_all(label, "\n")[[1]])
    mtext(label, at=mean(mp[i]), side=1, line = 0.5 + carriage_returns*0.5, cex=0.8)
  }

#   if (!is.null(plot_name))
#     text(max(mp)+1,95,plot_name, adj=c(1,1))
  dev.off()
}

data_out <- NULL
for (year in 1:length(datasets)) {
  year_dir <- file.path(temp_dir, names(datasets)[year])

  # now do the plot...
  data <- read.table(file.path(year_dir, "output.txt"), header=T)
  attribution_plot(file.path(output_dir, paste0(names(datasets)[year],".pdf")), data, col=col, plot_name=plot_names[year])
}

# combine plots
data_out <- NULL
for (year in 1:length(datasets)) {
  year_dir <- file.path(temp_dir, names(datasets)[year])

  data <- read.table(file.path(year_dir, "output.txt"), header=T)
  data$name <- plot_names[year]
  data_out <- rbind(data_out, data)
}


pdf(file.path(output_dir, "combined.pdf"), width=5, height=3.5)

d <- matrix(data_out$mean*100, nrow=length(datasets), byrow=T)
l <- matrix(data_out$lci*100, nrow=length(datasets), byrow=T)
u <- matrix(data_out$uci*100, nrow=length(datasets), byrow=T)
rownames(d) <- plot_names

ymax <- ceiling(max(d,u,l) / 20)*20
mp = barplot(d,beside=T, ylim=c(0,ymax),ylab="Percentage of cases attributed", cex.axis=0.8, cex.lab=0.8, legend=F)
segments(mp,l,mp,u,lwd=2)

for (i in 1:length(source_labels)) {
  label <- source_labels[as.character(i)]
  carriage_returns <- nrow(str_locate_all(label, "\n")[[1]])
  mtext(label, at=mean(mp[,i]), side=1, line = 0.5 + carriage_returns*0.5, cex=0.8)
}

dev.off()


