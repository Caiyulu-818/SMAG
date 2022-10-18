
# library
library(tidyverse)
library(rgdal)
library(ggplot2)
library(cowplot)
require(ggsci)
require(ggbeeswarm)


# files

mag.95 <- read_csv("dataset/mag95info.csv")
mag.all <- read_csv("/Users/caiyulu/Desktop/MAGcode/MAG/drep/mag40039_50_10.csv")
sample.info <- read_csv("dataset/sample_ecoinfo_r4.csv")

sample.info$ecosystem <- str_replace(sample.info$ecosystem, 
                                         pattern = "Agricultural L", 
                                         replacement = "Cultivated l") 

sample.info$ecosystem <- str_replace(sample.info$ecosystem, 
                                         pattern = "Water bodies",
                                         replacement = "Wetland")

# MAGs cluster size

quality1 <- rep("        36,393 medium quality MAGs\n", dim(mag.all)[1])
m<-quality1[mag.all$completeness > 90 & 
             mag.all$contamination < 5 & 
             mag.all$num_trna >=18 & 
             mag.all$`5S_5.8s` != 0 &
             mag.all$rna_16S!= 0 &
             mag.all$rrna_23S != 0
             ] <- 
    "\n 3646 high qulity MAGs"

pdf("/Users/caiyulu/Desktop/MAGcode/paper/Figure/fig1_pie.pdf",height = 3,width = 4)
par(mar = rep(0,4))
pie(table(quality1),cex = .7,
    radius = .5,
    init.angle = -105,
    col = pal_aaas(alpha = .75)(3)[2:1],
    border = FALSE)
dev.off()

gp.clu <- ggplot(mag.95, aes(x = cluster_members))+
        geom_point(stat = "count",shape = 1,size = 3, 
                   color = pal_aaas(alpha = .5)(3)[3]) +
        scale_x_log10(name = "Cluster size",
                      breaks = c(1,10,100),
                      labels = expression(10^0,10^1,10^2))+
        scale_y_log10(name = "Count",
                      breaks = c(1,10,100,1000,10000),
                      labels = expression(10^0,10^1,10^2,10^3,10^4))+
        theme_cowplot(font_size = 9)

gp.freq <- ggplot(sample.info, aes(x = Freq))+
        geom_point(stat = "count",shape = 1,size = 3, 
                   color = pal_aaas(alpha = .5)(3)[1]) +
        scale_x_log10(name = "MAGs numbers",
                      breaks = c(1,10,100,1000),limits = c(1,1000),
                      labels = expression(10^0,10^1,10^2,10^3))+
        scale_y_log10(name = "Count",
                      breaks = c(1,10,100,1000,10000),limits = c(1,1000),
                      labels = expression(10^0,10^1,10^2,10^3,10^4))+
        theme_cowplot(font_size = 9)

plot_grid(gp.clu, gp.freq,nrow = 2)
ggsave("figures/fig1_dis.pdf", height = 3,width = 2.5)
# make global map ######

# read shapefile
wmap <- readOGR(dsn="ne_110m_land", layer="ne_110m_land")
wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)
# convert to dataframe
wmap_df <- fortify(wmap)

# create a blank ggplot theme
theme_opts <-list(theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_rect(fill="#97CBFF"),
                        panel.border = element_blank(),
                        axis.line = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        plot.title = element_text(size=22,hjust = .5)))

grat <- readOGR("ne_110m_graticules_all", layer="ne_110m_graticules_15") 
grat_df <- fortify(grat)

bbox <- readOGR("ne_110m_graticules_all", layer="ne_110m_wgs84_bounding_box") 
bbox_df<- fortify(bbox)

grat_robin <- spTransform(grat, CRS("+proj=robin"))  # reproject graticule
grat_df_robin <- fortify(grat_robin)
bbox_robin <- spTransform(bbox, CRS("+proj=robin"))  # reproject bounding box
bbox_robin_df <- fortify(bbox_robin)

sample.info <- read_csv("sample_ecoinfo2.csv")
sel <- sample.info$lon > -120 & 
        sample.info$lon < -100 &
        sample.info$lat> 0 & 
        sample.info$lat < 30

sample.info <- sample.info[!sel,]

places_robin_df <- project(cbind(sample.info$lon[], 
                                 sample.info$lat[]),
                           proj="+init=ESRI:54030") 
places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")

places_robin_df <- places_robin_df %>% 
        mutate(ecosystem = sample.info$ecosystem,
               sample = sample.info$name,
               freq = sample.info$Freq) %>% 
        arrange(sample)

places_robin_df$ecosystem <- str_replace(places_robin_df$ecosystem, 
                                         pattern = "Agricultural L", 
                                         replacement = "Cultivated l") 

places_robin_df$ecosystem <- str_replace(places_robin_df$ecosystem, 
                                         pattern = "Water bodies",
                                         replacement = "Wetland")

sample.class <- sort(table(places_robin_df$ecosystem),decreasing = TRUE)
leg.lab <- paste0(names(sample.class)," (",sample.class,")")


gp.map <- ggplot(bbox_robin_df, aes(long,lat, group=group)) + 
        geom_polygon(fill=NA) +
        geom_polygon(data=wmap_df_robin, aes(long,lat, group=group, fill=hole)) + 
        #geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL),
        #          linetype="dashed", color="grey30",size =.1) +
        coord_equal() + 
        theme_opts + theme(plot.background = element_blank())+
        geom_point(data= na.omit(places_robin_df),
                   shape = 1,
                   aes(LONGITUDE, LATITUDE, 
                       color = ecosystem, size= freq),
                   inherit.aes = FALSE)+
        scale_size_continuous(range = c(0.3,8),name = "Number of MAGs", 
                              breaks = c(1,10,100,500),limits = c(1,600))+
        scale_color_npg(name = "Ecosystem\n(Number of samples)",
                                         breaks = names(sample.class), 
                                         labels = leg.lab)+
        scale_fill_manual(values=c("grey70", "#97CBFF40"), guide="none") +
        theme_void()+guides(colour = guide_legend(nrow = 3),alpha = "none")+
        theme(legend.position = "bottom",
              legend.text = element_text(size = 7),
              legend.title = element_text(size = 8, face = "bold"),
              legend.box = "vertical",legend.box.just = "left",
              legend.key.size = unit(1,"mm"),
              legend.spacing = unit(1,"mm"))

ggsave("figures/fig1_map.pdf",gp.map,height = 6,width = 7)


# MAGs quality

quality <- rep("Medium", dim(mag.all)[1])
quality[mag.all$completeness > 90 & 
            mag.all$contamination < 5 & 
            mag.all$num_trna >=18 & 
            mag.all$`5S(bac)/5.8s(arc)` != 0 &
            mag.all$`16S` != 0 &
            mag.all$`23S` != 0] <- "High"

## complete v.s. contamination

ggplot(mag.all, aes(y = contamination, x = completeness, color = quality))+
    geom_point()

min(mag.all$length[quality == "High"])

quality2 <- rep("Medium", dim(mag.all)[1])
quality2[mag.all$completeness > 90 & 
            mag.all$contamination < 5] <- "c"

a <- mag.all[mag.all$length == 532712 & quality == "High",]
b <- mag.all[mag.all$length == 12296985 & quality == "High",]


mag.entire <- mag.all %>% 
        select(length,completeness,contamination,N50,strain_heterogeneity) %>% 
        mutate(quality = quality)


g.len <- ggplot(mag.entire, aes(x = quality,
                                y = length,
                                color = quality)) +
        #geom_quasirandom(size = .01,
        #                 alpha = .01) +
        geom_quasirandom(alpha = .1, size = .5,
                         data = mag.entire[sample(1:23825,500),])+
        geom_boxplot(fill = NA,
                     color = "black",
                     outlier.size = 0) +
        theme_bw(base_size = 9) + coord_flip() +
        theme(
                aspect.ratio = .075,
                line = element_blank(),
                rect = element_blank(),
                axis.line.x = element_line(colour = "black"),
                axis.ticks.x = element_line(colour = "black"),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank()
        ) + ylab("Genome size (bp)") +guides(color = "none") +
        scale_y_log10(
                breaks = c(
                        seq(1e+5, 9e+5, by = 1e+5),
                        seq(1e+6, 9e+6, by = 1e+6),
                        seq(1e+7, 9e+7, by = 1e+7)
                ),
                limits = c(1e+5, 2e+7),
                labels = c(
                        expression(10 ^ 5),
                        rep("", 8),
                        expression(10 ^ 6),
                        rep("", 8),
                        expression(10 ^ 7),
                        rep("", 8)
                )
        )+
    scale_color_manual(values = pal_aaas()(2)[2:1])

g.complete <- ggplot(mag.entire[mag.entire$completeness >= 50 , ],
                     aes(x = quality,
                         y = completeness,
                         color = quality)) +
        #geom_quasirandom(size = .01,
        #                 alpha = .01) +
        geom_quasirandom(alpha = .1, size = .5,
                         data = mag.entire[sample(1:23825,500),])+
        geom_boxplot(fill = NA,
                     color = "black",
                     outlier.size = 0) +
        theme_bw(base_size = 9) + coord_flip() +guides(color = "none") +
        theme(
                aspect.ratio = .075,
                line = element_blank(),
                rect = element_blank(),
                axis.line.x = element_line(colour = "black"),
                axis.ticks.x = element_line(colour = "black"),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank()
        ) + ylab("Completeness (%)")+
    scale_color_manual(values = pal_aaas()(2)[2:1])


g.conta <- ggplot(mag.entire,
                  aes(x = quality,
                      y = contamination,
                      color = quality)) +
        #geom_quasirandom(size = .01,
        #                 alpha = .01) +
        geom_quasirandom(alpha = .1, size = .5,
                         data = mag.entire[sample(1:23825,500),])+
        geom_boxplot(fill = NA,
                     color = "black",
                     outlier.size = 0) +
        theme_bw(base_size = 9) + coord_flip() +guides(color = "none") +
        theme(
                aspect.ratio = .075,
                line = element_blank(),
                rect = element_blank(),
                axis.line.x = element_line(colour = "black"),
                axis.ticks.x = element_line(colour = "black"),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank()
        ) + ylab("Contamination (%)")+
        scale_color_manual(values = pal_aaas()(2)[2:1])

g.n50 <- ggplot(mag.entire,  aes(x = quality,
                                 y = N50,
                                 color = quality)) +
        #geom_quasirandom(size = .01,
        #                 alpha = .01) +
        geom_quasirandom(alpha = .1,size = .5, 
                         data = mag.entire[sample(1:23825,500),])+
        geom_boxplot(fill = NA,
                     color = "black",
                     outlier.size = 0) +
        theme_bw(base_size = 9) + coord_flip() +
        theme(
                aspect.ratio = .075,
                line = element_blank(),
                rect = element_blank(),
                axis.line.x = element_line(colour = "black"),
                axis.ticks.x = element_line(colour = "black"),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank()
        ) + ylab("N50 (bp)") +guides(color = "none") +
        scale_y_log10(
                breaks = c(
                        seq(1e+3, 9e+3, by = 1e+3),
                        seq(1e+4, 9e+4, by = 1e+4),
                        seq(1e+5, 9e+5, by = 1e+5),
                        seq(1e+6, 9e+6, by = 1e+6)
                ),
                limits = c(1e+3, 2e+6),
                labels = c(
                        expression(10 ^ 3),
                        rep("", 8),
                        expression(10 ^ 4),
                        rep("", 8),
                        expression(10 ^ 5),
                        rep("", 8),
                        expression(10 ^ 6),
                        rep("", 8)
                )
        )+
    scale_color_manual(values = pal_aaas()(2)[2:1])

g.hetero <- ggplot(mag.entire, aes(x = quality,
                                   y = strain_heterogeneity,
                                   color = quality)) +
        geom_quasirandom(data = mag.entire[sample(1:23825,500),] , 
                         inherit.aes = FALSE,size = .5,
                         alpha = .1,
                         aes(x = quality,
                             y = strain_heterogeneity,
                             color = quality))+
        geom_boxplot(fill = NA,
                     color = "black",
                     outlier.size = 0) +
        theme_bw(base_size = 9) + coord_flip() + guides(alpha = "none") +
        theme(
                aspect.ratio = .075,
                legend.position = "bottom",
                line = element_blank(),
                rect = element_blank(),
                legend.key.size = unit(2,"mm"),
                legend.text = element_text(size = 9),
                legend.text.align = 0,
                legend.title = element_text(size = 9,face = "bold"),
                #legend.direction = "vertical",
                axis.line.x = element_line(colour = "black"),
                axis.ticks.x = element_line(colour = "black"),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank()
        ) + ylab("Strain heterogeneity (%)")+
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))+
    scale_color_manual(values = pal_aaas()(2)[2:1],name =" MAG quality",
                         labels = c("High quality","Medium quality")
                       )

plot_grid(g.complete, g.conta, g.len, g.n50, g.hetero,
          rel_heights = c(1,1,1,1,2),
          align = "v", ncol = 1)

ggsave(filename = "figures/fig1_qual.pdf",height = 4, width = 3)

###### quality ######

mag.arc <- mag.all %>% 
    filter(phylum == "p__Patescibacteria") %>% 
    group_by(phylum,class) %>% 
    summarise(num = table(class)) %>% 
    arrange(desc(phylum))

mag.qual <- mag.all %>% 
    mutate(quality) %>% 
    select(length, contamination, completeness, quality) %>% 
    arrange(desc(quality))


qual.1 <-
    ggplot(mag.qual, aes(y = length, x = contamination)) +
    geom_point(alpha = .1, size = .5, aes(color = quality)) +
    coord_flip()+
    ylab("Genome size") + xlab("Contamination") +
    scale_y_log10(
        breaks = c(1e+5, 3e+5, 1e+6, 3e+6, 1e+7),
        limits = c(2.5e+5, 2e+7),
        labels = expression(1 * x * 10 ^ 5,
                            3 * x * 10 ^ 5,
                            1 * x * 10 ^ 6,
                            3 * x * 10 ^ 6,
                            1 * x * 10 ^ 7)
    ) +
    theme_cowplot(font_size = 9) +
    scale_color_manual(values = pal_aaas()(2)[2:1]) +
    theme(legend.position = "none", aspect.ratio = 1)

qual.2 <-
    ggplot(mag.qual, aes(x = length, y = completeness)) +
    geom_point(alpha = .1, size = .5,aes(color = quality)) + 
    xlab("Genome size") + ylab("Completeness") +
    scale_x_log10(
        breaks = c(1e+5, 3e+5, 1e+6, 3e+6, 1e+7),
        limits = c(2.5e+5, 2e+7),
        labels = expression(1 * x * 10 ^ 5,
                            3 * x * 10 ^ 5,
                            1 * x * 10 ^ 6,
                            3 * x * 10 ^ 6,
                            1 * x * 10 ^ 7)
    ) +
    theme_cowplot(font_size = 9) +
    scale_color_manual(values = pal_aaas()(2)[2:1]) +
    theme(legend.position = "none", aspect.ratio = 1)

qual.3 <-
    ggplot(mag.qual, aes(y = contamination, x = completeness, color = quality)) +
    geom_point(alpha = .1, size = .5) +
    ylab("contamination") + xlab("Completeness") +
    theme_cowplot(font_size = 9) +
    scale_color_manual(values = pal_aaas()(2)[2:1],name =" MAG quality",
                       labels = c("High quality","Medium quality")) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))+
    theme(legend.position = "bottom", aspect.ratio = 1,legend.key.size = unit(3,"mm"))

gp.qual <- plot_grid(
    qual.2,
    qual.1,
    qual.3,
    nrow = 3,
    rel_heights = c(1, 1, 1, 2),
    align = "hv",
    axis = "rlbt",
    labels = "auto"
)

ggsave("figures/qual.pdf",height = 7,width = 2.5)

mag.actino <- mag.95 %>% 
    filter(phylum == "p__Actinobacteriota") %>% 
    filter(class == "c__")

