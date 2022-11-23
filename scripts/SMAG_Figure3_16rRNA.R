# MAG 16S rRNA

rrna1 <- read_table("dataset/magotu_tax.txt",col_names = F)

rrna2 <- rrna1 %>% 
        filter(str_detect(X2, pattern = "Actinobacteria")) %>% 
        mutate(phyla = str_sub(X2, 
                               start = str_locate(X2,"p:")[,2]+1, 
                               end = str_locate(X2,"c:")[,1]-1)) %>% 
        mutate(similar = str_sub(phyla, 16,21)) %>% 
        arrange(similar)

rrna.tree <- read.tree(file = "dataset/mafft_actin16s_bmge.fasta.treefile")
rrna.tree$node.label <- 1:4836
rrna.tree.class <- rrna.tree$tip.label
novel <- rep("a",4838)
novel[str_detect(rrna.tree.class, "_tax_")] <- "b" 

class <- str_sub(rrna.tree.class,
                 start = str_locate(rrna.tree.class,"_c_")[,2]+1,
                 end = str_locate(rrna.tree.class,"_o_")[,1]-1)

x <- tibble(label = rrna.tree$tip.label, trait = novel, class)
rrna.tree1 <- full_join(rrna.tree, x, by="label")

clade <- tree_subset(rrna.tree1, node=7438, levels_back=0)

rrna.p <- ggtree(clade, layout = "rect",branch.length = "none") 


rrna.p %>% 
        scaleClade(265, 20) %>% 
        scaleClade(140, .5) %>% 
        collapse(140, "max", fill = pal_aaas(alpha = .75)(7)[3]) %>% 
        collapse(222, "max", fill = pal_aaas(alpha = .75)(7)[5]) %>% 
        collapse(211, "max", fill = pal_aaas(alpha = .75)(7)[4]) %>% 
        collapse(241, "max", fill = pal_aaas(alpha = .75)(7)[6]) %>% 
        collapse(253, "max", fill = pal_aaas(alpha = .75)(7)[2]) %>% 
        collapse(262, "max", fill = pal_aaas(alpha = .75)(7)[7]) %>% 
        collapse(265, "max", fill = pal_aaas(alpha = .75)(7)[1]) 

ggsave("figures/rrna.pdf", ,width = 2,height = 4)
