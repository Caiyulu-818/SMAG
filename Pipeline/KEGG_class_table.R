## making KEGG class table

kegg.json <- read_tsv("dataset/ko00001.keg",skip = 3, col_names = FALSE)

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

kegg.pathway.table <- tibble(class.A.col1, class.B.col1, class.C.col1, ko, class.D.value1)

write.csv(x = kegg.pathway.table, "kegg_level_table.csv")