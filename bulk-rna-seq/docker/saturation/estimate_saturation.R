#!/usr/bin/env Rscript


library('preseqR')
library('reshape2')
library('ggplot2')
library('optparse')

option_list <- list(
  make_option(c('-i','--input'), type="character", default=NULL,
              help="input file name", metavar="character"),
  make_option(c('-n','--name'), type="character", default=NULL,
              help="library name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)

opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop('Input file not specified')
}

# input params
input_file <- opt$input
lib_name <- opt$name

#lib_name <- 'saturation_curve'
#input_file <- 'dedup_stats_per_umi_per_position.tsv'

d <- read.csv(input_file, sep='\t')
d2 <- as.matrix(d[,c(1,2)])

# Get the total library size
lib_size = sum(d$counts * d$instances_pre)

# Sequence of values to calculate expected molecules for
s <- seq(0.5,10,0.1)

m <- data.frame(
  count = lib_size * s,
  preseqR.rSAC = preseqR.rSAC(d2)(s)
)

write.csv(m,'saturation.tsv')

m2 <- melt(m,id.vars = c('count'))

actual_sequencing_point = data.frame(
  count = lib_size * 1,
  preseqR.rSAC = preseqR.rSAC(d2)(1)
)

p <- ggplot(m2, aes(x=count, y=value, color=variable)) + 
  geom_line() + geom_abline(intercept = 0, slope=1) + 
  geom_point(data = actual_sequencing_point, mapping=aes(x=count, y = preseqR.rSAC), inherit.aes = FALSE) +
  scale_x_continuous(limits=c(0,max(m$count))) + 
  scale_y_continuous(limits=c(0,max(m$count))) +
  ggtitle(paste0(lib_name, ' Saturation Curve')) +
  theme_bw() + theme(legend.position = "none")

ggsave(filename = "saturation_curve.pdf", plot = p, width=7,height=7)
