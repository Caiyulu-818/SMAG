# dataset loading

snv.count <- read_table("snv_dataset/snv_pos_count.tsv", col_names = c("b","a")) %>% 
        arrange(desc(b)) %>% 
        mutate(cumm = cumsum(b)) %>% 
        mutate(sample = 1:2358)

### cummulated 
ggplot(snv.count, aes(x = sample, y = cumm)) +
        xlab("Number of species")+ylab("Cummulated number of SNVs")+
        theme_cowplot(font_size = 9)+
        geom_hline(yintercept = 180538439, linetype = 2)+
        geom_line(size = 1, color = pal_aaas()(2)[2])+
        scale_y_log10(labels = expression(1.0*x*10^7,1.0*x*10^8,1.8*x*10^8),breaks = c(1e+7,1e+8,1.8e+8))+
        theme(aspect.ratio = .618)


# dataset
#mag.all <- read_csv("dataset/mag40039_50_10.csv")

# SNV selection
ani.tab <- fread("dataset/Cdb.csv")
snv.sampleinfo <- readxl::read_xlsx("snv_dataset/snvcatalog_maginfo.xlsx",sheet = 1,col_names = c("a","b")) 

snv.sampleinfo1 <- snv.sampleinfo %>% 
        filter(!is.na(a)) %>% 
        mutate(size = mag.all$length[pmatch(b, mag.all$user_genome)]) %>%
        mutate(phylum = mag.all$phylum[pmatch(b, mag.all$user_genome)]) %>% 
        mutate(hetero = mag.all$strain_heterogeneity[pmatch(b, mag.all$user_genome)]) %>% 
        mutate(complete = mag.all$completeness[pmatch(b, mag.all$user_genome)]) %>% 
        group_by(a) %>% 
        mutate(gsize = sum(100*size/complete)) %>% 
        filter(!duplicated(a)) %>% 
        mutate(cluster = ani.tab$secondary_cluster[pmatch(paste0(b,".fa") , table = ani.tab$genome,duplicates.ok = T)]) %>% 
        mutate(novel = mag.tree$novel[pmatch(str_remove(cluster,"_[0-9]+$"), str_remove(mag.tree$cluster, pattern = "_[0-9]+$"),duplicates.ok = T)])

snv.mean <- as.tibble(merge(snv.count, snv.sampleinfo1, by = "a")) %>% 
        filter(!is.na(novel)) %>% 
        group_by(a) %>% 
        mutate(gsize = mean(size)) %>% 
        mutate(snv.density = b.x/gsize*1000) %>%
        ungroup() %>% 
        group_by(phylum, novel) %>% 
        mutate(pmean = median(snv.density,na.rm = T)) %>%  
        mutate(pmean1 = median(b.x, na.rm = T)) %>%  
        ungroup() %>% 
        arrange(desc(b.x)) %>% 
        mutate(cumm =  cumsum(b.x)) %>% 
        group_by(novel) %>% 
        mutate(cumm1 = cumsum(b.x))

snv.text <- tibble(x = 2100,y=c(1.3e+8,4.7e+7,7.6e+7),label = c("SGBs","kSGBs","uSGBs"))

ggplot(na.omit(snv.mean), aes(x = sample, y = cumm1, color = novel, group = novel)) +
        xlab("Number of species")+ylab("Cummulated SNV numbers")+
        theme_cowplot(font_size = 10)+
        geom_hline(yintercept = c(143652279,57148685,86503594), linetype = 2)+
        geom_line(size = 1) +
        geom_line(size = 1, aes(x = sample, y = cumm),color = pal_aaas()(3)[3],inherit.aes = F) +
        scale_color_aaas()+
        geom_text(aes(x,y,label = label),data = snv.text,inherit.aes = F,size = 3)+
        scale_y_continuous(labels = expression(0,1.0*x*10^7,1.0*x*10^8,1.4*x*10^8,5.7*x*10^7,8.6*x*10^7),
                           breaks = c(0,1e+7,1e+8,1.4e+8,5.7e+7,8.6e+7))+
        theme(legend.position = "none")

ggsave("figures/snv_curve.pdf",height = 3, width = 4)


#write.csv(snv.mean, "dataset/snv_dentisy.csv")

ggplot(na.omit(snv.mean), aes(
        x = fct_reorder(fct_lump(str_remove(phylum,"p__"), 30), pmean),
        y = snv.density,
        color = novel
                )) +
        #geom_boxplot(outlier.shape = 1, color = pal_aaas(alpha = .75)(5)[5]) + #geom_violin(fill = NA)+
        stat_summary() +
        scale_y_log10(breaks = c(1e-3,1e-2,1e-1,1e+0,1e+1,1e+2),name = "# SNVs per kb",
                      labels = expression(10^-3,10^-2,10^-1,10^0,10^1,10^2)) +
        coord_flip() + xlab("")+
        scale_color_aaas(alpha = .5)+
        theme_cowplot(font_size = 9) +
        theme(legend.position = "none")

ggsave("figures/snv_phyla.pdf",height = 5, width = 3)

snv.mean$novel <- factor(snv.mean$novel, labels = c("kSGB","uSGB"))

ggplot(na.omit(snv.mean), aes(x = gsize, y = snv.density)) +
        geom_point(size = 2, shape = 1, alpha = .5) +
        scale_y_log10(breaks = c(1e-3,1e-2,1e-1,1e+0,1e+1,1e+2),name = "# SNVs per kb",
                      labels = expression(10^-3,10^-2,10^-1,10^0,10^1,10^2)) +
        scale_x_log10(limits = c(3e+5,1e+7),name = "Genome size (Mb)",
                      breaks = c(3e+5,1e+6,3e+6,1e+7),
                      labels = c("0.3","1.0","3.0","10.0")) +
        theme_cowplot(font_size = 9) +
        facet_grid(~ novel) +
        geom_smooth(method = "lm", se = FALSE, color = "darkred",size = 1) +
        theme(legend.position = "none",panel.spacing = unit(3, "mm"))

ggsave("figures/snv_gsize.pdf",height = 2.5, width = 3.5)



ggplot(na.omit(snv.mean), aes(x = hetero, y = snv.density)) +
        geom_point(size = 2, shape = 1, alpha = .5) +
        scale_y_log10(breaks = c(1e-3,1e-2,1e-1,1e+0,1e+1,1e+2),name = "# SNVs per kb",
                      labels = expression(10^-3,10^-2,10^-1,10^0,10^1,10^2)) +
        theme_cowplot(font_size = 9) +
        facet_grid(~ novel) + xlab("Stain heterogeneitiy")+
        geom_smooth(method = "lm", se = FALSE, color = "darkred",size = 1) +
        theme(legend.position = "none",aspect.ratio = 1)

ggsave("figures/snv_hetero.pdf",height = 2, width = 4)


ggplot(na.omit(snv.mean), aes(x = gsize, y = snv.density)) +
        geom_point(size = 2, shape = 1, alpha = .5) +
        scale_y_log10(breaks = c(1e-3,1e-2,1e-1,1e+0,1e+1,1e+2),name = "# SNVs per kb",
                      labels = expression(10^-3,10^-2,10^-1,10^0,10^1,10^2)) +
        scale_x_log10(limits = c(3e+5,1e+7),name = "Genome size (Mb)",
                      breaks = c(3e+5,1e+6,3e+6,1e+7),
                      labels = c("0.3","1.0","3.0","10.0")) +
        theme_cowplot(font_size = 9) +
        facet_wrap(~ fct_lump(str_remove(phylum, "p__"), 19), nrow = 5) +
        geom_smooth(method = "lm", se = FALSE, color = "darkred",size = 1) +
        theme(legend.position = "none",aspect.ratio = .5,panel.spacing = unit(3, "mm"))

ggsave("figures/snv_phylo_seq.pdf",height = 5, width = 6)


