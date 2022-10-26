# functional annotation
## packages
require(tidyverse)
require(KEGGREST)
require(ggsci)
require(cowplot)
require(jsonlite)

## making KEGG class table

kegg.json <- read_tsv("/Users/caiyulu/Desktop/MAGcode/function/kegg_csv/ko00001.keg",skip = 3, col_names = FALSE)

kegg.df <- kegg.json %>% 
        filter(X1 != "B") %>% 
          mutate(class = str_sub(X1,1,1))

kegg.ad <- kegg.df %>% 
        filter(class %in% c("A","D"))

class.A <- grep("A",kegg.ad$class)
class.A.value <- kegg.ad$X1[class.A]
class.A.col <- rep(class.A.value, (c(class.A[-1], 58118) - class.A)-1)
class.A.col1 <- str_sub(class.A.col, start = 8, end = str_length(class.A.col))

kegg.bd <- kegg.df %>% 
        filter(class %in% c("B","D"))

class.B <- grep("B  ",kegg.bd$X1)
class.B.value <- kegg.bd$X1[class.B]
class.B.col <- rep(class.B.value, (c(class.B[-1], 58166)-class.B)-1)
class.B.col1 <- str_sub(class.B.col, start = 10, end = str_length(class.B.col))


kegg.cd <- kegg.df %>% 
        filter(class %in% c("C","D"))

class.C <- grep("C",kegg.cd$class)
class.C.value <- kegg.cd$X1[class.C]
class.C.col <- rep(class.C.value, (c(class.C[-1], 58647)-class.C)-1)
class.C.col1 <- str_sub(class.C.col, start = 12, end = str_length(class.C.col)-15)

class.D <- grep("D",kegg.df$class)
class.D.value <- kegg.df$X1[class.D]
ko <- str_sub(class.D.value, 8, 13)
class.D.value1 <- str_sub(class.D.value, 16, str_length(class.D.value))
gene.name <- str_sub(class.D.value1, 1, 
                     str_locate(class.D.value1,"\\;")[,1] - 1)

ec.number <- str_sub(class.D.value1, 
                     str_locate(class.D.value1,"EC:")[,1], 
                     str_length(class.D.value1))

kegg.pathway.table <- tibble(class.A.col1, class.B.col1, class.C.col1, ko, class.D.value1,gene.name,ec.number)

write.csv(x = kegg.pathway.table, "kegg_level_table.csv")

# annotation results  #########
kegg.phylo <- read_csv("dataset/KEGG_annotation_MAGs.csv")
mag.all <- read_csv("dataset/mag40039_50_10.csv")

kegg.phyla.df <- kegg.phylo %>%
        mutate(complete = mag.all$completeness[pmatch(user_genome, mag.all$user_genome,duplicates.ok = TRUE)]) %>% 
        mutate(contamin = mag.all$contamination[pmatch(user_genome, mag.all$user_genome,duplicates.ok = TRUE)]) %>% 
        mutate(length = mag.all$length[pmatch(user_genome, mag.all$user_genome,duplicates.ok = TRUE)])%>% 
        mutate(class.b = kegg.pathway.table$class.B.col1[pmatch(function_type,
                                                                kegg.pathway.table$class.C.col1,
                                                                duplicates.ok = T)]) %>% 
        filter(complete >= 90 & contamin <= 5) 

kegg.phyla.df1 <- kegg.phyla.df %>% 
        group_by(class.b, phylum, info) %>%
        filter(!is.na(info)) %>%
        summarise(pathway.count = sum(count)) %>%
        spread(info, pathway.count) %>%
        ungroup() %>% 
        mutate(class.a = kegg.pathway.table$class.A.col1[pmatch(class.b,
                                                                kegg.pathway.table$class.B.col1,
                                                                duplicates.ok = T)]) 

gene.occur <- rep("Both", 1446)
gene.occur[is.na(kegg.phyla.df1$known)] <- "uSGB"
gene.occur[is.na(kegg.phyla.df1$unknown)] <- "kSGB"

kegg.phyla.df2 <- kegg.phyla.df1 %>% 
        mutate(gene.occur = gene.occur) %>% 
        filter(!class.a %in% c("Human Diseases","Organismal Systems")) %>% 
        filter(!is.na(class.a)) %>% 
        group_by(gene.occur) %>% 
        mutate(nu = n())

ggplot( kegg.phyla.df2,
        aes(    x = fct_infreq(class.b),
                y = fct_infreq(str_remove(phylum,pattern = "p__")),
                #fill = gene.occur,
                color = gene.occur,#alpha = nu,
                shape = gene.occur
        )) +
        geom_point(size = 3) +
        scale_color_aaas(name = "") +
        #scale_fill_aaas(name = "") +
        scale_shape_manual(name = "",values = c(22,22,15))+
        theme_cowplot(font_size = 9) +
        ylab("Phylum") + 
        xlab("Functional class") +
        coord_flip()+
        facet_grid(class.a~., scales = "free", space = "free") +
        theme(strip.text.y = element_text(angle = 0,face = "bold"),
              axis.title = element_text(face = "bold"),
              axis.text.x = element_text(angle = 45,
                                         hjust = 1),
              legend.position = "top")

 ggsave("figures/kegg_phyla.pdf",height = 3.5, width = 14)

# COG

cog.df <- read_csv("dataset/cog_all_sub_un.csv")

cog.df1 <- cog.df %>% 
        group_by(COG_category, info) %>% 
        summarise(count = sum(sum)) %>% 
        ungroup()

require(treemapify)

treemapify(cog.df1, area = "count",subgroup = "COG_category")

ggplot(na.omit(cog.df1), aes(area =  count, 
                    fill = info, 
                    subgroup = COG_category))+
        geom_treemap()+
        geom_treemap_text(aes(label = COG_category))+
        geom_treemap_subgroup_border()


# pangenomes ###############

pangenome.mag <- readxl::read_xlsx("dataset/pan_anno_yjw/pairwise_maginfo.xlsx", 
                                   col_names = c("x1","x2"))

unknown.info <-  read_csv("snv_dataset/known_maginfo.csv", 
                          col_names = c("num","V1","ani","af","novel")) %>% 
        arrange("V1") %>% 
        mutate(genome = str_remove(V1, ".report.genbank")) %>% 
        mutate(genome = str_remove(genome, ".report.other")) %>% 
        mutate(genome = str_remove(genome, ".report")) %>% 
        mutate(genome = str_remove(genome, ".fa")) 

pangenome.mag.df <- pangenome.mag %>% 
        mutate(phylum = mag.all$phylum[pmatch(str_remove(pangenome.mag$x2,"\\."), 
                                              str_remove(mag.all$user_genome,"\\."),
                                              duplicates.ok = TRUE)]) %>% 
        mutate(length = mag.all$length[pmatch(str_remove(pangenome.mag$x2,"\\."), 
                                              str_remove(mag.all$user_genome,"\\."),
                                              duplicates.ok = TRUE)]) %>% 
        mutate(novel = unknown.info$novel[pmatch(str_remove(pangenome.mag$x2,"\\."), 
                                              str_remove(unknown.info$genome,"\\."),
                                              duplicates.ok = TRUE)]) %>% 
        group_by(x1) %>% 
        mutate(gnum = n()) %>% 
        mutate(len1 = max(length,na.rm = T)) %>% 
        filter(!is.na(novel)) %>% 
        filter(!duplicated(x1))

file.names <- list.files("dataset/pan_anno_yjw/list/")

gene.total.number <- NULL

for(i in file.names){
        gf <- dim(read_csv(paste0("dataset/pan_anno_yjw/list/",i)))[1]
        gene.total.number <- c(gene.total.number, gf)
}

total.gene.number <- tibble(pgenome = str_sub(file.names, 1, 5),
                            gene = str_sub(file.names, 7, 10),
                            gnum = gene.total.number)

core.gene.number <- total.gene.number %>% 
        filter(gene == "core"& gnum != 0)

acc.gene.number <- total.gene.number %>% 
        filter(gene == "acce"& gnum != 0)

pan.gene.df <- merge.data.frame(core.gene.number, acc.gene.number, by = "pgenome") %>% 
        as_tibble() %>% 
        mutate(percent = gnum.x/(gnum.y+gnum.x)*100) 

pan.gene.df1 <- merge(pan.gene.df, pangenome.mag.df, by.x = "pgenome",by.y = "x1")

gp1 <- ggplot(na.omit(pan.gene.df1), aes(y = percent, x = gnum))+
        geom_point(shape = 1, size =3,aes(color = phylum))+
        coord_cartesian(clip = "off") +
        ylab("% core gene")+ xlab("# of genomes")+ #xlim(0,55)+
        theme_cowplot(font_size = 10)+
        theme(legend.position = "none")+
        geom_smooth(method = "glm", se = F,color = "darkred")

gp2 <- ggplot(na.omit(pan.gene.df1), aes(y = percent, x = len1))+
        geom_point(shape = 1, size =3,aes(color = phylum))+
        #coord_cartesian(clip = "off") +
        ylab("% core gene")+ xlab("Genome size(Mbp)")+ #xlim(0,55)+
        theme_cowplot(font_size = 10)+
        theme(legend.position = "none")+
        scale_x_log10(labels = c(1,3,10),
                      breaks = c(1e+6,3e+6,1e+7))+
        geom_smooth(method = "lm", se = F,color = "darkred")

plot_grid(gp1, gp2,align = "h",nrow = 1,labels = "auto")
ggsave("figures/pangenome_core_size_genome.pdf",height = 3,width = 6)

require(ggridges)

ggplot(na.omit(pan.gene.df1), 
       aes(x = fct_reorder(str_remove(phylum,"p__"), percent), 
           y = percent)) +
        geom_point(aes(color = novel, shape = novel),
                   size = 4)+
        coord_flip()+ xlab("")+ ylab("% core genes")+
        scale_shape_manual(values = c(1, 16), name = "",
                           labels = c("kSGB", "uSGB"))+
        scale_color_manual(name = "", 
                           values = pal_aaas(alpha = .75)(3)[2:3],
                           labels = c("kSGB", "uSGB"))+
        theme_cowplot(font_size = 10)+
        theme(legend.position = c(.6,.2))

ggsave("figures/pangenome_core_proportion.pdf",
       height = 4,width = 3)

ggplot(na.omit(pan.gene.df1), 
       aes(y = fct_reorder(str_remove(phylum,"p__"), percent))) +
        geom_density_ridges2(aes(fill = phylum, x = percent), 
                            panel_scaling = FALSE,
                            stat = "binline",
                            alpha = .5,
                            color = "white",
                             scale = 2) +
        #scale_y_discrete(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0))  +
        coord_cartesian(clip = "off") +
        theme_ridges(grid = FALSE,font_size = 9) + 
        xlab("Proportion of core gene (%)")+ 
        ylab("")+ xlim(0,100)+
        theme(legend.position = "none")



### core and accessory genes annotated with different database

core.go <- read_csv("dataset/result_v2/Core.gene2go.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "core") %>% 
        mutate(db = "GO") %>% 
        mutate(per = gnum/core.gene.number$gnum[pmatch(X1,table = core.gene.number$pgenome)])

acc.go <- read_csv("dataset/result_v2/Accessory.gene2go.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "acc") %>% 
        mutate(db = "GO") %>% 
        mutate(per = gnum/acc.gene.number$gnum[pmatch(X1,table = acc.gene.number$pgenome)])


core.ko <- read_csv("dataset/result_v2/Core.gene2ko.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "core") %>% 
        mutate(db = "KEGG") %>% 
        mutate(per = gnum/core.gene.number$gnum[pmatch(X1,table = core.gene.number$pgenome)])


acc.ko <- read_csv("dataset/result_v2/Accessory.gene2ko.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "acc") %>% 
        mutate(db = "KEGG") %>% 
        mutate(per = gnum/acc.gene.number$gnum[pmatch(X1,table = acc.gene.number$pgenome)])


core.cazy <- read_csv("dataset/result_v2/Core.gene2CAZy.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "core") %>% 
        mutate(db = "CAZy") %>% 
        mutate(per = gnum/core.gene.number$gnum[pmatch(X1,table = core.gene.number$pgenome)])

acc.cazy <- read_csv("dataset/result_v2/Accessory.gene2CAZy.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "acc") %>% 
        mutate(db = "CAZy") %>% 
        mutate(per = gnum/acc.gene.number$gnum[pmatch(X1,table = acc.gene.number$pgenome)])

core.cog <- read_csv("dataset/result_v2/Core.gene2COG.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "core") %>% 
        mutate(db = "COG") %>% 
        mutate(per = gnum/core.gene.number$gnum[pmatch(X1,table = core.gene.number$pgenome)])


acc.cog <- read_csv("dataset/result_v2/Accessory.gene2COG.csv", skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "acc") %>% 
        mutate(db = "COG") %>% 
        mutate(per = gnum/acc.gene.number$gnum[pmatch(X1,table = acc.gene.number$pgenome)])

core.egg <- read_csv("dataset/result_v2/Core.gene2eggnog.csv",skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "core") %>% 
        mutate(db = "eggNOG") %>% 
        mutate(per = gnum/core.gene.number$gnum[pmatch(X1, table = core.gene.number$pgenome)])

acc.egg <- read_csv("dataset/result_v2/Accessory.gene2eggnog.csv",skip = 1,col_names = F) %>% 
        group_by(X1, X2) %>% 
        summarise(gnum = n()) %>%
        group_by(X1) %>% 
        summarise(gnum = n()) %>%
        mutate(gene = "acc") %>% 
        mutate(db = "eggNOG") %>% 
        mutate(per = gnum/acc.gene.number$gnum[pmatch(X1,table = acc.gene.number$pgenome)])

gene.num.df <- add_row(core.egg, acc.egg) %>%
        add_row(core.go) %>%
        add_row(core.ko) %>%
        add_row(core.cazy) %>%
        add_row(core.cog) %>%
        add_row(acc.cog) %>%
        add_row(acc.ko) %>%
        add_row(acc.go) %>%
        add_row(acc.cazy) %>%
        mutate(novel = pangenome.mag.df$novel[pmatch(X1, pangenome.mag.df$x1, duplicates.ok = TRUE)])




wilcox.fun <- function(x = gene.num.df){
        x1 = x[x$novel == "known", 5]
        x2 = x[x$novel == "unknown", 5]
        wil.p <- t.test(as.numeric(x1$per), as.numeric(x2$per),alternative = "greater")$p.value
        return(wil.p)
}

wilcox.fun(gene.num.df[gene.num.df$gene == "core" & gene.num.df$db == "eggNOG", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "core" & gene.num.df$db == "COG", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "core" & gene.num.df$db == "KEGG", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "core" & gene.num.df$db == "GO", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "core" & gene.num.df$db == "CAZy", ])

wilcox.fun(gene.num.df[gene.num.df$gene == "acc" & gene.num.df$db == "eggNOG", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "acc" & gene.num.df$db == "COG", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "acc" & gene.num.df$db == "KEGG", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "acc" & gene.num.df$db == "GO", ])
wilcox.fun(gene.num.df[gene.num.df$gene == "acc" & gene.num.df$db == "CAZy", ])

gene.num.pval <- gene.num.df %>% 
        group_by(gene,db) %>% 
        summarise(pval = wilcox.fun(.data))


gene.num.df$gene <- factor(gene.num.df$gene, levels = c("core","acc"),
                           labels = c("Accessory genes","Core genes")[2:1])
ggplot(gene.num.df, aes(x = fct_relevel(db,"eggNOG","COG","KEGG","GO","CAZy"), 
                        y = per*100, 
                        fill = novel))+
        #stat_summary()+
        geom_violin(draw_quantiles = .5,scale = "width",width = .75, size = .5)+
        #scale_y_log10(limit = c(0.01,100),name= "% annotated genes",
        #              breaks = c(.01,.1,1,10,100),
        #              labels = expression(10^-2,10^-1,10^0,10^01,10^02))+  
        xlab("")+ ylim(0,100)+ ylab("Annotated\nproporiton(%)")+
        scale_fill_manual(breaks = c("known","unknown"),
                         values = pal_aaas(alpha = .75)(3)[2:3],
                         name = "",
                         labels = c("kSGB","uSGB"))+
                         #labels = c("Core genes","Accessary genes")) +
        #geom_quasirandom(size = 1, alpha =.5, shape = 1)+
        facet_grid(~gene)+
        theme_cowplot(font_size = 10)+
        theme(legend.position = c(.85,.9))

ggsave("figures/pangenome_gene_annotation.pdf",height = 2,width = 7)

core.ko.table <- read_csv("dataset/result_v2/Core.gene2ko.csv", skip = 1,col_names = F) %>% 
        mutate(ko = str_sub(X3, 4, 9)) %>% 
        mutate(class.c = kegg.pathway.table$class.B.col1[pmatch(ko , kegg.pathway.table$ko, duplicates.ok = T)]) %>% 
        group_by(X1, class.c) %>% 
        summarise(path.num = n()) %>% 
        mutate(gene = 1) %>% 
        ungroup()

acc.ko.table <- read_csv("dataset/result_v2/Accessory.gene2ko.csv", 
                         skip = 1,col_names = F) %>% 
        mutate(ko = str_sub(X3, 4, 9)) %>% 
        mutate(class.c = kegg.pathway.table$class.B.col1[pmatch(ko , kegg.pathway.table$ko, duplicates.ok = T)]) %>% 
        group_by(X1, class.c) %>% 
        summarise(path.num = n())%>% 
        mutate(gene = 2) %>% 
        ungroup()

ko.table <- add_row(acc.ko.table, core.ko.table) %>% 
        spread(key = class.c, value = path.num, fill = 0) %>% 
        mutate(novel = pangenome.mag.df$novel[pmatch(X1, pangenome.mag.df$x1, duplicates.ok = TRUE)])

ko.table.usgb <- ko.table %>% 
        filter(novel == "unknown") %>% 
        group_by(X1) %>% 
        mutate(num = n()) %>% 
        filter(num == 2)

ko.table.ksgb <- ko.table %>% 
        filter(novel == "known") %>% 
        group_by(X1) %>% 
        mutate(num = n()) %>% 
        filter(num == 2)

cohen.value.usgb <- tibble()
for(i in c(3:10,12:15, 18:33,35:42,44:45,47:55)){
        d1 <- psych::cohen.d(ko.table.usgb[, c(2,i)], group = "gene")
        cohen.value.usgb <- rbind(cohen.value.usgb ,d1$cohen.d)
}

cohen.value.ksgb <- tibble()
for(i in c(3:10,12:15, 18:33,35:42,44:45,47:55)){
        d1 <- psych::cohen.d(ko.table.ksgb[, c(2,i)], group = "gene")
        cohen.value.ksgb <- rbind(cohen.value.ksgb ,d1$cohen.d)
}

cohen.value.usgb <- tibble()
for(i in c(3:10,12:15, 18:33,35:42,44:45,47:55)){
        d1 <- psych::cohen.d(ko.table.usgb[, c(2,i)], group = "gene")
        cohen.value.usgb <- rbind(cohen.value.usgb ,d1$cohen.d)
}

wilcox.value.ksgb <- NULL
for(i in c(3:10,12:15, 18:33,35:42,44:45,47:55)){
        wil <- t.test(unlist(ko.table.ksgb[seq(1,74,2), c(i)]), 
                           unlist(ko.table.ksgb[seq(2,74,2), c(i)]),
                           alternative = "t",paired = TRUE)
        wilcox.value.ksgb <- rbind(wilcox.value.ksgb ,wil$p.value)
}

wilcox.value.usgb <- NULL
for(i in c(3:10,12:15, 18:33,35:42,44:45,47:55)){
        wil <- wilcox.test(unlist(ko.table.usgb[seq(1,74,2), c(i)]), 
                           unlist(ko.table.usgb[seq(2,74,2), c(i)]),
                           alternative = "t",paired = TRUE)
        wilcox.value.usgb <- rbind(wilcox.value.usgb ,wil$p.value)
}


cohen.df1 <- as_tibble(cohen.value.ksgb) %>% 
        mutate(pval = wilcox.value.ksgb[,1]) %>% 
        mutate(pathway = row.names(cohen.value.ksgb)) %>% 
        mutate(class = kegg.pathway.table$class.A.col1[pmatch(pathway ,kegg.pathway.table$class.B.col1)]) %>% 
        filter(!class %in% c("Human Diseases", "Organismal Systems", "Brite Hierarchies"))

cohen.df2 <- as_tibble(cohen.value.usgb) %>% 
        mutate(pval = wilcox.value.usgb[,1]) %>% 
        mutate(pathway = rownames(cohen.value.usgb)) %>% 
        mutate(class = kegg.pathway.table$class.A.col1[pmatch(pathway , kegg.pathway.table$class.B.col1)]) %>% 
        filter(!class %in% c("Human Diseases", "Organismal Systems", "Brite Hierarchies"))

cohen.df <- add_row(cohen.df1, cohen.df2) %>% 
        mutate(novel = c(rep("kSGB",29), rep("uSGB",29))) %>% 
        #filter(pval < 0.1) %>% 
        arrange(effect, desc(novel)) %>% 
        filter(pathway != "") %>% 
        filter(pval < 0.05) %>% 
        filter(!is.na(class))

cohen.df3 <- cohen.df %>% 
        filter(novel == "kSGB") %>% 
        arrange(effect)

ggplot(cohen.df, 
       aes(x = fct_reorder(pathway, effect), 
                     y = effect, 
           fill = class)) +
        geom_bar(stat = "identity") +
        coord_flip(ylim = c(0,2)) + 
        ylab("Effect size") + xlab("")+
        theme_cowplot(font_size = 10)+
        facet_grid(~novel, scale = "free",space = "free") +
        scale_y_continuous(breaks = c(0,1,2)) +
        theme(legend.position = "right")+
        scale_fill_manual(name = "", 
                          guide = guide_legend(nrow = 5,
                                               keyheight = unit(1.5,"mm")),
                          values = c(pal_aaas(alpha = .75)(4),"grey"))


ggsave("figures/pangenome_enrich_pathway.pdf", height = 3, width = 7)



ggplot(cohen.df, 
       aes(x = fct_relevel(pathway, unique(pathway)), 
           color = novel)) +
        geom_segment(aes(y = lower,yend = upper,
                         xend = fct_relevel(pathway, unique(pathway)))) +
        geom_point(aes(x = fct_relevel(pathway, unique(pathway)),
                       y = effect))+
        coord_flip()+ ylab("Effect size")+xlab("")+
        theme_cowplot(font_size = 10)+
        theme(legend.position = c(.6,.2),strip.text.y = element_text(angle = 0))+
        scale_color_manual(name = "", values = pal_aaas(alpha = .75)(3)[2:3])

# Actinobacteria function


actino.kegg <- read_csv(unzip("dataset/p_actin_ko.csv.zip"), 
                        col_names = FALSE, skip = 1)

nit.kegg <- read_csv("dataset/actin_ntri_ko.csv",col_names = FALSE, skip = 1) %>% 
        mutate(class.c = kegg.pathway.table$class.C.col1[pmatch(X3,
                                                        kegg.pathway.table$ko,
                                                        duplicates.ok = T)]) %>% 
        mutate(class.a = kegg.pathway.table$class.A.col1[pmatch(X3,
                                                                kegg.pathway.table$ko,
                                                                duplicates.ok = T)]) %>%
        mutate(class.b = kegg.pathway.table$class.B.col1[pmatch(X3,
                                                                kegg.pathway.table$ko,
                                                                duplicates.ok = T)]) %>% 
        filter(class.a %in% c("Metabolism")) %>%
        group_by(class.c, class.b) %>% 
        summarise(sum = n()) %>% 
        mutate(class = "Nitriliruptoria") 

#ko.table <- read_csv("dataset/keggko.csv")

ko.actino1 <- actino.kegg %>%
        mutate(class.c = kegg.pathway.table$class.C.col1[pmatch(X3,
                                                                kegg.pathway.table$ko,
                                                                duplicates.ok = T)]) %>% 
        mutate(class.a = kegg.pathway.table$class.A.col1[pmatch(X3,
                                                                kegg.pathway.table$ko,
                                                                duplicates.ok = T)]) %>%
        mutate(class.b = kegg.pathway.table$class.B.col1[pmatch(X3,
                                                                kegg.pathway.table$ko,
                                                                duplicates.ok = T)]) %>% 
        filter(class.a %in% c("Metabolism")) %>%
        mutate(class = mag.all$class[pmatch(X2,
                                            mag.all$user_genome,
                                            duplicates.ok = T)]) %>% 
        group_by(class.c, class.b, class) %>% 
        summarise(sum = n()) 

ko.actino2 <- ko.actino1 %>% 
        ungroup() %>% 
        add_row(nit.kegg) %>%
        filter(!str_detect(class, "UBA"))%>%
        filter(!str_detect(class, "RBG")) %>% 
        filter(class.c != "G") %>% 
        mutate(class = str_replace(class, "c__$",
                                   "Candidate Gramensolibacteria")) %>% 
        mutate(class = str_remove(class, "c__"))

ggplot(ko.actino2, aes(y = fct_infreq(class.c), 
                       x = fct_infreq(class))) +
        geom_tile(aes(fill = class.c, color = class),size = .25) +
        scale_color_aaas()+
        theme_classic(base_size = 10) +
        facet_grid(class.b ~ ., scale = "free",space = "free") +
        theme(axis.text.y = element_text(size = 4),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(angle = 45, 
                                         size = 10, 
                                         hjust = 1),
              strip.text.y = element_text(angle = 0),
              legend.position = "none")+xlab("")+ylab("")

ggsave("figures/Actinobacteriota_function.pdf", height = 10,width = 6)

# CRISPR ####

crispr.info <- read_table("dataset/crispr2mag.txt",col_names = FALSE) %>% 
        mutate(bin.name = str_sub(X1, start = 1, end = str_locate(crispr.info$X1, "_")[,1]-1)) %>% 
        group_by(bin.name) %>% 
        summarise(freq = n())

