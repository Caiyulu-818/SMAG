# library
require(ggtree)
require(ggtreeExtra)
require(tidytree)
require(tidyverse)
require(treeio)
library(ggstar)
library(ggnewscale)
library(treemap)
library(colorspace)
library(ggsci)
library(cowplot)

#sample.info <- read_csv("sample_ecoinfo2.csv")

unknown.info <-  read_csv("/Users/caiyulu/Desktop/MAGcode/MAG/unknown/ref_mag/known_maginfo.csv") %>% 
        arrange("V1")

cultured.info <-  read_csv("/Users/caiyulu/Desktop/MAGcode/MAG/unknown/ref_mag/cultured_info_refani95.csv") %>% 
        arrange(MAG)
mag.95<-read.csv("/Users/caiyulu/Desktop/MAGcode/MAG/drep/mag95info.csv")
bac<-subset(mag.95,classification=="d__Bacteria")

tree.file <- read.tree("/Users/caiyulu/Desktop/MAGcode/paper/Figure/bact.tre")
tree.file<-read.tree("/Users/caiyulu/Desktop/MAGcode/MAG/new_genome/tree/RAxML_bestTree.mag20177prodigal_refined1.tre")
       # read.tree("dataset/RAxML_bestTree.mag20177prodigal_refined1.tre")

sample.info<-read.csv("/Users/caiyulu/Desktop/MAGcode/paper/data_submit/submit/Supplemental information/supplymental_tables_r1/TableS1_related_toFigure1.xlsx")

mag.tree <- mag.95 %>%
        arrange(genome) %>% 
        mutate(eco = sample.info$ecosystem[pmatch(
                str_remove(mag.95$genome,pattern = "r*bin.[0-9]*.fa"),
                sample.info$mag_name,
                duplicates.ok = TRUE
        )]) %>%
        mutate(phy = phylum == "p__") %>% 
        mutate(culture = cultured.info$status) %>% 
        mutate(novel = unknown.info$...5) %>% 
        filter(genome %in% tree.file$tip.label) 

usgb <- rep("z",19187)
usgb[mag.tree$species == "s__"] <- "s"
usgb[mag.tree$genus == "g__"] <- "g"
usgb[mag.tree$family == "f__"] <- "f"
usgb[mag.tree$order == "o__"] <- "o"
usgb[mag.tree$class == "c__"] <- "c"
usgb[mag.tree$phylum == "p__"] <- "p"

mag.tree <- mag.tree %>% 
    mutate(usgb)

tree.df <- tree.file  %>% as_tibble() %>%
        full_join(mag.tree, by = c("label" = "genome")) %>%
        filter(!is.na(node)) %>%
        as.treedata()

edge.colors <- rep(alpha("black", alpha = .25), 39720)
edge.colors[tree.df@data$novel == "unknown"] <-  pal_aaas(alpha = .75)(5)[3]
edge.colors[tree.df@data$culture == "cultured"] <-  pal_aaas()(5)[2]


tree.df2 <- tibble(
        node = 1:39720,
        color = edge.colors
)

set.seed(1234)
phy.color = qualitative_hcl(20)[sample(1:20,replace = F)]

gptree <- ggtree(
        tree.df,
        layout = "fan",
        ladderize = TRUE,
        open.angle = 20,
        size = .05
        ) %<+% tree.df2 +  
        aes(color = I(color)) +
    new_scale_color()+new_scale_fill()+
    geom_fruit(
        geom = geom_tile,
        pwidth = 0.2,
        width = 0.06,
        offset = 0,
        #axis.params = list(title = "Domain",size = 4),
        mapping = aes(
            y = tip,
            fill = classification,
            color = classification
        )) +
    scale_fill_aaas() +
    scale_color_aaas() +
    new_scale_color()+new_scale_fill()+
        geom_fruit(
                geom = geom_tile,
                pwidth = 0.2,
                width = 0.06,
                offset = 0.02,
                mapping = aes(
                        y = tip,
                        fill = fct_lump(str_remove(phylum,"p__"),19),
                        color = fct_lump(str_remove(phylum,"p__"),19)
        )) +
    scale_fill_manual(values = phy.color) +
    scale_color_manual(values = phy.color) +
    new_scale_color()+new_scale_fill()+
        geom_fruit(
                geom = geom_tile,
                pwidth = .2,
                width = 0.1,
                offset = 0.03,
                mapping = aes(
                        y = tip,
                        x = fct_relevel(usgb, "s","g","f","o","c","p","z"),
                        color = usgb,
                        fill = usgb
                )
        ) +
        scale_color_manual(values = c(pal_aaas()(6),"white")) +
        scale_fill_manual(values = c(pal_aaas()(6),"white"))  +
        new_scale_color()+new_scale_fill()+
        geom_fruit(
                geom = geom_bar,
                orientation = "y",
                offset = .03,
                stat = "identity",
                width = .01,
                color = pal_aaas(alpha = .5)(5)[5],
                fill = pal_aaas(alpha = .5)(5)[5],
                mapping = aes(y = tip,
                              x = size)
        ) +
  new_scale_color()+new_scale_fill()+
        geom_fruit(
                geom = geom_bar,
                orientation = "y",
                stat = "identity",
                offset = .01,
                width = .01,
                color = pal_aaas(alpha = .5)(5)[4],
                fill = pal_aaas(alpha = .5)(5)[4],
                mapping = aes(y = tip,
                              x = log10(cluster_members))
        )+
    theme(legend.key.height = unit(2,units = "mm"))

ggsave(gptree,
       filename = "/Users/caiyulu/Desktop/MAGcode/paper/Figure/fig2_tree1.pdf",
       height = 10,
       width = 10)

# make archaea and bacteria seperate

arch.tip <- mag.tree$genome[mag.tree$classification == "d__Archaea"]
bact.tip <- mag.tree$genome[mag.tree$classification != "d__Archaea"]

arch.tree <- drop.tip(tree.df, bact.tip) %>% as.phylo()
bact.tree <- drop.tip(tree.df, arch.tip) %>% as.phylo()

write.tree(arch.tree, file = "/Users/caiyulu/Desktop/MAGcode/paper/Figure/arch.tre")
write.tree(bact.tree, file = "/Users/caiyulu/Desktop/MAGcode/paper/Figure/bact.tre")

# make rarefy curve

## phylum
rare.phy <- NULL
for(i in seq(1,40000,1000)){
    for(j in 1:10){
        rare.phy <- c(rare.phy,length(unique(sample(mag.all$phylum,i))))
    }
    print(i)
}

rare.class <- NULL
for(i in seq(1,40000,1000)){
    for(j in 1:10){
        rare.class <- c(rare.class,length(unique(sample(mag.all$class,i))))
    }
    print(i)
}

rare.order <- NULL
for(i in seq(1,40000,1000)){
    for(j in 1:10){
        rare.order <- c(rare.order,length(unique(sample(mag.all$order,i))))
    }
    print(i)
}

rare.family <- NULL
for(i in seq(1,40000,1000)){
    for(j in 1:10){
        rare.family <- c(rare.family,length(unique(sample(mag.all$family,i))))
    }
    print(i)
}

rare.genus <- NULL
for(i in seq(1,40000,1000)){
    for(j in 1:10){
        rare.genus <- c(rare.genus,length(unique(sample(mag.all$genus,i))))
    }
    print(i)
}

rare.species <- rarefy(mag.95$cluster_members,sample = seq(1,40000,1000),se = F)

rare.phy.df <- data_frame(i = c(rep(rep(seq(1,40000,1000),each = 10),5),seq(1,40000,1000)),
                          rare = c(rare.phy,rare.class,rare.order,rare.family,rare.genus,c(rare.species)),
                          taxa = c(rep(c("Phylum","Class","Order","Family","Genus"),each = 400),rep("Speices",40)))  %>% 
    group_by(i,taxa) %>% 
    summarise(r = mean(rare))

rare.phy.df$taxa <- factor(rare.phy.df$taxa, levels = c("Phylum","Class","Order","Family","Genus","Speices"))


ggplot(rare.phy.df, aes(x = i, y = r, color = taxa)) +
    geom_line(size = .2) +
    theme_classic()+xlab("Random selected MAGs")+ylab("Taxa number")+
    scale_color_aaas(name = "")+theme(legend.key.height = unit(2,"mm"))

ggsave("figures/fig2_rare.pdf",height = 2,width = 4)

## uSGB and kSGB in phylum
mag.tree$novel[mag.tree$culture == "cultured"] <- "cultured"

sgb.df1 <- mag.tree %>% 
        select(novel, phylum, classification,cluster_members) %>% 
        group_by(novel,classification, phylum) %>% 
        summarise(clu = sum(cluster_members)) %>% 
        mutate(phy = str_remove(phylum,"p__")) %>% 
        mutate(dom = str_remove(classification,"d__")) %>% 
        ungroup() %>% 
        group_by(phy) %>% 
        mutate(sum = sum(clu)) 


gp.sgb1 <- ggplot(sgb.df1[sgb.df1$sum > 270,], 
                  aes(x = fct_reorder(phy, sum), 
                       y= clu, 
                       fill = fct_relevel(novel, c("cultured","known","unknown")[3:1])))+
               geom_bar(stat = "identity") +
        facet_grid(fct_relevel(dom,c("Bacteria","Archaea"))~.,
                   scales = "free",space = "free")+
        coord_flip()+
        ylab("Number of MAGs") + 
        xlab("")+ 
        theme_cowplot(font_size = 9)+
        scale_fill_manual(name = "", labels = c("Cultivated kSGB","Uncultivated kSGB","uSGB")[3:1],values = pal_aaas()(3)[3:1])+
        theme(legend.position = c(.3, .3),
              axis.text.y = element_text(face = "italic"),
              strip.placement = "outside",strip.text.y = element_text(angle = 0))

ggsave("/Users/caiyulu/Desktop/MAGcode/paper/Figure/fig2_SGBs.pdf",height = 3,width = 4)

# uSGB and kSGB in various samples

sgb.df2 <- mag.tree %>% 
        select(novel, eco, cluster_members) %>% 
        group_by(novel,eco) %>% 
        summarise(clu = sum(cluster_members)) %>% 
        ungroup() %>% 
        group_by(eco) %>% 
        mutate(sum = sum(clu))

gp.sgb2 <- ggplot(na.omit(sgb.df2), aes(x = fct_reorder(eco,sum), 
                                            y= clu, 
                                           fill = fct_relevel(novel,c("unknown","known","cultured"))))+
        geom_bar(stat = "identity")+
        coord_flip()+
        ylab("Number of MAGs")+ 
        xlab("")+
        theme_cowplot(font_size = 9)+
        scale_fill_aaas(name = "", labels = c("uSGB","Uncultivated kSGB","Cultivated kSGB"))+
        theme(legend.position = c(.5, .5),
              axis.text.y = element_text(face = "italic"),
              strip.placement = "outside")

plot_grid(gp.sgb2, gp.sgb1,rel_heights = c(1,3),nrow = 2,
          align = 'hv',axis = "tblr",label_x = .08,
          labels = c("b","c"))
ggsave("figures/fig2_SGBs.pdf")

seq.dis <- sample.info[!is.na(sample.info$Freq), ] %>% 
    arrange(ecosystem) %>% 
    filter(Freq < 1000)

ggplot(seq.dis,
       aes(x = number_seqs, y = Freq,color = ecosystem))+
        geom_point(alpha = .1, shape =16,size = 2)+
        scale_x_log10(labels = expression(10^6,10^7,10^8,10^9),
                      breaks = c(1e+6,1e+7,1e+8,1e+9),
                      #limits = c(1e+6,1e+9),
                      name = "Number of reads after quaility control") +
        geom_smooth(method = "lm",se = F,color = 'darkred',size = .5,linetype = 2)+
        scale_y_log10(limits = c(1,400),
                      name = "Reconstructed MAGs per sample") +
        facet_wrap(~fct_infreq(ecosystem),scale = "free_y") +
        theme_cowplot(font_size = 9) +
        theme(aspect.ratio = .5,legend.position = "none")
ggsave(filename = "figures/seqdep.pdf",height = 4,width = 6)


clust.df <- tibble(mag = c("14072 singleton SBGs (70.9%)",
                           "2777 SGBs with 2 MAGs (14.0%)",
                           #"1254 SGBs with 3 MAGs(6.3%)",
                           "3012 SGBs > 3 MAGs (15.1%)"),
                   num = c(14072,2777,3012)) #1254,1758))

require(ggpubr)
pdf("figures/pie.pdf",height = 2,width = 8)
par(mar=c(0,8,0,8))
pie(clust.df$num,
    labels = clust.df$mag, 
    cex = 1,
    init.angle = 50,
    radius = 1,font =2,
    col = pal_aaas(alpha = .8)(4),
    border = NA)
dev.off()

ggplot(mag.tree,aes(x = fct_reorder(fct_lump(phylum, 19),size), y = size,
                    fill = novel)) +
    geom_violin(draw_quantiles = .5)+
    scale_y_log10() +
    theme_cowplot(font_size = 10)+
    coord_flip()+
    scale_color_discrete_divergingx()+
    theme(legend.position = "")


mag.cluster <- mag.tree %>% 
    group_by(cluster_members,novel) %>% 
    summarise(freq = n())

ggplot(mag.cluster, aes(x = cluster_members, y = freq, fill =  novel)) +
    geom_bar(stat = "identity")+
    scale_fill_aaas()




