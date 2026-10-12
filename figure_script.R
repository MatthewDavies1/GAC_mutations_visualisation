
library(dplyr)
library(ggplot2)
library(gggenes)
library(grid)
library(lubridate)

raw.data <- read.csv("annotated_snps_input_R.csv", header = TRUE, sep = ";")

input.data <- raw.data[,c(1,7,2,3,5,11,13,18,19)]

genes <- data.frame(
  operon = rep("gac operon", 12),
  genename = c("gacA","gacB","gacC","gacD","gacE","gacF","gacG","gacH","gacI","gacJ","gacK","gacL"),
  start = c(604904,605877,607021,607956,608759,609989,610993,612735,615388,616085,616419,617686),
  end = c(605758,607031,607953,608759,609964,610996,612738,615209,616083,616426,617705,619182)
)
operonarrows <- ggplot(genes, aes(xmin = start, xmax = end, y = operon, fill = genename, label = genename)) +
  geom_gene_arrow(arrow_body_height = unit(8, "mm"),
                  arrowhead_height = unit(8, "mm"),
                  arrowhead_width = unit(2, "mm")) + 
  geom_gene_label(height = unit(4, "mm"), fontface = "italic") +
  xlim(c(604766,619182)) +
  theme_genes() + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none")
operonarrows

#synonymous mutations
input.data.synonymous <- input.data %>% 
  filter(grepl("synonymous_variant", S.NS_1) | grepl("synonymous_variant", S.NS_2)) %>%
  mutate(sum = (ifelse(S.NS_1 == "synonymous_variant", freq, 0) + ifelse(S.NS_2 == "synonymous_variant", freq_2,0)))
input.data.synonymous <- input.data.synonymous %>%
  mutate(relabund = (sum/2021)*100)

synonymous <- ggplot(input.data.synonymous, aes(x = pos, y = relabund)) +
  geom_col(width = 5, alpha = 0.9, position = "identity") +
  xlim(c(604766,619182)) + 
  scale_y_continuous(labels=function(x) paste0(x,"%"), limits = c(0,100)) +
  theme_classic() +
  ylab("0") +
  ggtitle("Synonymous")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18, colour = "white"),
        axis.text.y = element_text(size = 12))
synonymous

#nonsynonymous mutations
input.data.nonsyn <- input.data %>% 
  filter(grepl("missense_variant", S.NS_1) | grepl("missense_variant", S.NS_2) | grepl("stop", S.NS_1)) %>%
  mutate(sum = (ifelse(S.NS_1 == "missense_variant", freq, 0) + ifelse(S.NS_2 == "missense_variant", freq_2,0)))
input.data.nonsyn <- input.data.nonsyn %>%
  mutate(relabund = (sum/2021)*100)

nonsynonymous <- ggplot(input.data.nonsyn, aes(pos, relabund)) +
  geom_col(width = 5, alpha = 0.9, position = "identity") +
  xlim(c(604766,619182)) + 
  scale_y_continuous(labels=function(x) paste0(x,"%"), limits = c(0,100)) +
  theme_classic() + 
  ylab("SNP frequency (n = 2021)")+
  ggtitle("Nonsynonymous")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 12))
nonsynonymous

#intergenic mutations
input.data.intergenic <- input.data %>% 
  filter(grepl("upstream", S.NS_1) | grepl("upstream", S.NS_2)) %>%
  mutate(sum = (ifelse(S.NS_1 == "upstream_gene_variant", freq, 0) + ifelse(S.NS_2 == "upstream_gene_variant", freq_2,0)))
input.data.intergenic <- input.data.intergenic %>%
  mutate(relabund = (sum/2021)*100)

intergenic <- ggplot(input.data.intergenic, aes(pos, relabund)) +
  geom_col(width = 5, alpha = 0.9, position = "identity") +
  xlim(c(604766,619182)) + 
  scale_y_continuous(labels=function(x) paste0(x,"%"), limits = c(0,100)) +
  theme_classic() + 
  ylab("0")+
  ggtitle("Intergenic")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18, colour="white"),
        axis.text.y = element_text(size = 12))
intergenic

grid.newpage()
snpfigure <-grid.arrange(synonymous,nonsynonymous,intergenic,operonarrows, ncol = 1, heights = c(1,1,1,0.2))
ggsave("SNP figure for illustrator- FINAL.png", snpfigure, dpi = 500, width = 14, height = 10, units = "in")
