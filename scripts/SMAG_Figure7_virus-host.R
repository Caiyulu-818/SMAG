# viral-host
require(tidyverse)
require(igraph)
require(ggraph)
require(tidygraph)
require(ggstar)

## loading files
virus.file <- read_csv("dataset/viruslinkages_MAG.csv")
mag.95 <- read_csv("dataset/mag95info.csv")
mag.all <- read_csv("dataset/mag40039_50_10.csv") 

unknown.info <-  read_csv("snv_dataset/known_maginfo.csv") %>% 
        arrange("V1") %>% 
        mutate(genome = str_remove(V1, ".report.genbank")) %>% 
        mutate(genome = str_remove(genome, ".report.other")) %>% 
        mutate(genome = str_remove(genome, ".report")) %>% 
        mutate(genome = str_remove(genome, ".fa")) 

virus.assign.family <- read_csv("dataset/family_pro.csv")
virus.assign.genus <- read_csv("dataset/genus_pro.csv")

cluster.info <- read_csv(file = "dataset/Cdb.csv")

mag.cluster <- cluster.info$primary_cluster[pmatch(mag.all$user_genome,
                                                   paste0(cluster.info$genome, ".fa"),
                                     duplicates.ok = TRUE)]

cluster.rep <- mag.95$genome[pmatch(mag.cluster, 
                                    str_remove(mag.95$cluster, "_[0-9]*&"),
                                    duplicates.ok = TRUE)]

mag.novel <- unknown.info$...5[pmatch(cluster.rep , 
                                      paste0(unknown.info$genome, ".fa"), 
                                      duplicates.ok = TRUE)] 

mag.novel[is.na(mag.novel)] <- "unknown"

mag.all <- mag.all %>% 
        mutate(novel = mag.novel)

# virus host 
virus.host.df <- virus.file %>%
        mutate(phylum = mag.all$phylum[pmatch(bin, 
                                              mag.all$user_genome, 
                                              duplicates.ok = T)]) %>%
        mutate(novel = mag.all$novel[pmatch(bin, 
                                            mag.all$user_genome, 
                                            duplicates.ok = T)])

virus.class.df <- virus.host.df %>% 
        group_by(phylum, novel, method) %>% 
        summarise(count = n()) %>% 
        group_by(phylum) %>% 
        mutate(sum = sum(count)) %>% 
        mutate(per = count/sum) %>% 
        filter(sum > 100)

ggplot(virus.class.df[!duplicated(virus.class.df$phylum), ], 
        aes(x = fct_reorder(str_remove(phylum,"p__"), sum), 
                           y = sum)) +
        geom_bar(stat = "identity", fill = pal_npg()(3)[3]) +
        geom_text(aes(label = sum),hjust = 0, size = 3,fontface = "bold")+
        scale_y_log10(breaks = c(1,10,100,1000,10000),
                      limits = c(1, 100000), 
                      labels = expression(10^0,10^1,10^2,10^3,10^4)) +
        ylab("Count") + xlab("Phylum") +
        coord_flip() +
        theme_classic()+
        theme(aspect.ratio = 1.5)

ggsave("figures/viral_abund.pdf",height = 3,width = 3)

ggplot(virus.class.df, aes(x = fct_reorder(str_remove(phylum, "p__"), sum), 
                           y = per*100, 
                           fill = novel)) +
        geom_bar(stat = "identity",position = "stack") +
        scale_fill_manual(values = pal_aaas()(3)[2:3],
                          name = "",
                          labels = c("kSGB", "uSGB"))+
        ylab("Percent (%)") + xlab("Phylum") +
        coord_flip() +
        theme_classic(base_size = 10)+
        theme(legend.position = "bottom",aspect.ratio = 1.5,
              legend.key.height = unit(2,"mm"))
ggsave("figures/viral_novel.pdf",height = 3,width = 3)

ggplot(virus.class.df, aes(x = fct_reorder(str_remove(phylum, "p__"), sum), 
                           y = per*100, 
                           fill = method)) +
        geom_bar(stat = "identity",position = "stack") +
        scale_fill_manual(name = "",values = pal_aaas()(5)[4:5],
                       labels = c("Spacer", "Prophage"))+
        ylab("Percent (%)") + xlab("Phylum") +
        coord_flip() +
        theme_classic(base_size = 10)+
        theme(legend.position = "bottom",aspect.ratio = 1.5,
              legend.key.height = unit(2,"mm"))
ggsave("figures/viral_method.pdf",height = 3,width = 3)

virus.host.df1 <- virus.host.df %>% 
        #filter(!duplicated(paste(virus,phylum))) %>% 
        group_by(virus) %>% 
        mutate(v.num = n()) %>% 
        filter(v.num > 25) 
      
virus.host.df1 <- na.omit(virus.host.df1)
  
ggplot(na.omit(virus.host.df1), 
               aes(x = fct_reorder(virus, phylum), 
                   y = fct_infreq(bin)))+
        geom_point(aes(color = novel)) +
        #geom_tile(aes(color = virus, 
        #              fill = virus))+
        theme_classic()+
        facet_grid(phylum~., scales = "free", space = "free")+
        guides(color ="none",
               fill = "none")

el <- na.omit(as.matrix(virus.host.df1[,c(2,4)]))

write.csv(el, "edge_list.csv")

g.vh <- graph_from_edgelist(el = el, directed = TRUE)
write_graph(g.vh1, file = "g_vh.graphml",format = "graphml")

node.phylum <- rep("virus", 1569)
node.phylum[pmatch(virus.host.df1$bin , V(g.vh)$name,duplicates.ok = T)] <- virus.host.df1$phylum
edge.novel <- rep("unknown", 2878)
edge.novel[pmatch(virus.host.df1$bin, el[,1], duplicates.ok = T)] <- virus.host.df1$novel


virus.host.df1[virus.host.df1$virus == "GSV_42450", ]
node.name <- V(g.vh)$name
node.name[!str_detect(node.name,pattern = "GSV")] <- ""


g.vh1 <- g.vh %>% 
        as_tbl_graph() %>% 
        activate("nodes") %>% 
        mutate(phylum = node.phylum) %>% 
        mutate(degree = degree(g.vh)) %>% 
        mutate(vname = node.name) %>% 
        activate("edges") %>% 
        mutate(novel = edge.novel) %>% 
        arrange(desc(novel))

node.type <- node.phylum
node.type[!str_detect(node.phylum,pattern = "virus")] <- "bacteria"

g.vh2 <- read.graph("gene_taxa1.graphml",format = "graphml") %>% 
        as_tbl_graph() %>% 
        activate("nodes") %>% 
        mutate(node.type )

as_tbl_graph(b) %>% 
        activate("nodes") %>% 
        filter(name %in% c[[1]]$name)

d <- degree(b)

d[c("GSV_39462", "GSV_66726","GSV_270")]

neighborhood(a, order = 1, nodes = "GSV_42450") -> c

mycolors <- colorRampPalette(pal_aaas()(9))(18)

ggraph(g.vh2, layout = "manual", x = x, y = y)+
        geom_edge_arc(alpha = .5, 
                      width = .25,
                      strength = .2,
                           aes(color = novel)) +
        scale_edge_color_manual(values = pal_aaas()(3)[2:3], 
                                name = "Host type",
                                labels = c("kSGB","uSGB")) +
        #geom_node_text(aes(label = vname),hjust = 1,size = 2)+
        new_scale_color() +
        geom_node_point(#shape = 16,
                        aes(size = degree, shape = node.type,
                            color = fct_lump(gsub("p__","",phylum), 18))) +
        scale_size_continuous(range = c(1,4))+
        scale_color_manual(values = mycolors) +
        theme_graph() +
        theme(aspect.ratio = 1)

ggraph(g.vh2) +
        geom_edge_diagonal(alpha = .25, width = .2,
                           aes(color = novel)) +
        scale_edge_color_manual(values = pal_aaas()(3)[2:3], 
                                name = "Host type",
                                labels = c("kSGB","uSGB")) +
        #geom_node_text(aes(label = vname),hjust = 1,size = 2)+
        new_scale_color()+
        geom_node_point(shape = 1,alpha = .5,
                        aes(size = degree, 
                        color = fct_lump(gsub("p__","",phylum), 9))) +
        scale_color_npg(name = "Host phylum")+
        scale_size_continuous(range = c(.5,5),breaks = c(10,50,100), name = "Degree")+
        coord_flip()+
        theme_graph(base_family = "sans",base_size = 10)+
        theme(legend.key.size = unit(2,"mm"))

ggsave("figures/virus_host_net.pdf", height = 8, width = 4)

