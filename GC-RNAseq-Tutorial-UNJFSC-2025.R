#### TUTORIAL DE RNA-seq usando edgeR
#### Gianmarco Castillo Huaccho
#### Codigo de Orcid: https://orcid.org/0000-0003-0865-5717
#### Github: https://github.com/GianmarcoCastillo
#### Biologo con mencion Biotecnologo egresado de UNJFSC
#### Este tutorial es Free, puede usarlo quien quiera sin olvidar citar el trabajo
#### El tutorial fue realizado con la version de R.4.4.3 - 15 de abril 2025 -

# Libraries
library(edgeR)
library(statmod)
library(ComplexHeatmap)
library(ggplot2)
y <- read.table("/home/marcos/Descargas/RNA-seq/RNAseq/drought.cq01.counts", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
head(y)

#write.csv(counts,"drought_cq01.csv")
#Defaults
### PRE control pre floraci?n
### D_PRE Control sequia 
####
#### PASO 2  Crear tratamientos y colores
####
treat<-factor(c(rep("PRE", 6), rep("C_POS", 3), rep("D_POS", 3)))
cols<-factor(c(rep("darkorange", 6), rep("purple", 3), rep("darkgreen", 3)))
#### 
#### PASO 3 Crear el objeto DGEList
####
cq01 <- DGEList(y, group = treat)
cq01$samples
### Variable de que genes me quedo y que genes no me quedo 
### la funcion row SUms a las filas le calcula la suma  
###  >=3 se selecuiona si cumple o no la condicion de que un gen se repita en mas de 3 muestras biologicas
#### LOs falsos no tienen mas de 20 lecturas en mas de 3 librerias. 
####
#### PASO 4 FILTRAR
####

keep <- rowSums(y>=20)>=3
tail(keep)
summary(keep)
f01 <- cq01 [keep, , keep.lib.size=FALSE]
f01$samples
f01

##### Si el numero normalizado de read es constante en una libreria normalizada
png(file="barplot.png",width = 25, height = 25, units = "cm", res = 300)
barplot(f01$samples$lib.size,names=colnames(f01$counts),las=2)

dev.off()
#### Normalizaremos las lecturas por 1 millon 
#### 0.01% = 100 en un millon; 10%= 100 
#### CPM CONTEOS A UN MILLON; LLUEGO SE SACA EL LOGARITMO BASE 10
boxplot(cpm(f01$counts, log = T),names=colnames(f01$counts),las=2)


#### pARA VER LAS MULTIDIMENCIONAL; Y nos ayuda a obtener en dos dimenciones 
### En este caso por mas que el color verde este separado hay la posibilidad que puedan tener patrones
png(file="MDS.png",width = 25, height = 25, units = "cm", res = 300)
plotMDS(f01, pch = c(rep(1, 6),rep(2, 3), rep(3, 3)), col=as.character(cols),cex=1.5, main="M-values")
dev.off()

####
#### PASO  6 NORMALIACION
####
#### Factor de normalizacion aplica favtores de normalizaci?n; del tercer cuartil 75%

N01 <- calcNormFactors(f01, robust=T) ### el robust es de 2 pasos

boxplot(cpm(N01, log = T),las=2)
plotMDS(N01, pch = c(rep(1, 6),rep(2, 3), rep(3, 3)), col=as.character(cols),cex=1.5, main="M-values")

plotMDS(N01, col=as.numeric(N01$samples$group))
#legend("bottomleft", as.character(unique(N01$samples$group)), col=1:3, pch=20)
####
#### PASO 7 EstimaciOn de la disperciOn
####

treat
### como los grupos; aqui puedes color las otras variables que influyen en el modelo
####  Factores de interferencia representa un tipo de organo. 
#model.matrix(~jaula+treat+tiempo)
design <- model.matrix(~0+treat)
#### 
design
#### Aqui se obtienen los valores para estimar la disperci?n de todos los valores que necesitamos para ajustar a
##la funciOn de lo que 
N01<- estimateDisp(N01, design = design, robust = T)
N01
### d1 <- estimateCommonDisp(N01, verbose=T )
### names(d1)
####
#### PASO 8 regresi?n 
####
### En este paso se calcula cuando tiene poco numeros de muestras
fit01 <- glmQLFit(N01, robust = T)
fit01
## Tratamiento uno lo trata como positivo, el tratamiento 2 lo trata como negativo (-1) y al tratamiento
### 3 no lo considere (0), siempr eienen que tener signos opuestos
#### Se comparo  $comparison "1*treatC_POS -1*treatD_POS"; Aqui puedes contrast = c(1,-1,0)) contrastar respecto a tus tratamientos
### log (T1/T2)= 2=2^1
### log (T2/T1)= 1/2= 2^-1
### logCPM Conteos por millon (lllevas todos tus read a un millon)
#### logFC logaritmo de Fold Change; cuantas veces esta repetido 
dge01 <- glmQLFTest(fit01, contrast = c(1,-1,0))
dge01

tags01 <- topTags(dge01, n=Inf) #### n = Inf es para seleccione a todos tus valores
#### FDR (falsos positivos; correccion del Pvalue usando bonferroni)es de 0.000594 (hay diferencias entre tratamientos )la probabilidad de este AUR62008337
####  2^4 = -4.230434 = y quiere decir que esta 16 veces uno respecto al otro tratamiento
#### De la muestra el p-value (cual es la probabilidad de obtener un gen al azar de una muestra; si es muy baja se encuentra en la cola de una distribuci?n de datos )
#### El criterio del umbral 
head(tags01$table)
tags01$comparison

### FDR <= 0.01
### logFC >= 2
###  lfc 2 es estatico si no que depende de lo que quieres buscar (gretty pondria 4 por que quiere encontrar y optimizar la identificacion de genes celulasa)

sig <- decideTests(dge01, p.value = 0.01, lfc = 2)
summary(sig)
# write.table(tags01$table, file = "CPOSvsDPOS.DGE.tsv", sep="\t", quote = F, col.names = T, row.names = T)
sel <- tags01$table[abs(tags01$table$logFC) >= 1 & tags01$table$FDR <=  0.01, ]
dim(sel)

#sel <- abs(tags01$table$logFC) > 1 &tags01$table$FDR <=  0.01
#summary(sel)

#head(tags01$table)
##### Esto es para hacer True osea son genes elejidos 
sel <- tags01$table$FDR <= 0.01 & abs (tags01$table$logFC) >= 1
head(sel)

##### en este paso el decide test 
sig<- decideTests(dge01, p.value = 0.01, lfc=1, p.adjust="BH")
head(sig)
summary(sig)
####
####  Obtencion de la tabla normalizadas de los datos de los genes diferenciamente expresados
####
logcounts <- cpm (N01, log=T)
head(logcounts)
genes <- rownames(sig)
s.g. <- genes [sig<0 | sig>0]
### decimos que seleccione a los genes GDE entre las 
s.lcpm <- logcounts[s.g., 7:12]
head(s.lcpm)
dim(s.lcpm)
Heatmap(as.matrix(s.lcpm))

#######################
####################### PLOT HEATMAP
#######################
#######################
logcounts <- cpm (N01, log=T)
sig<- decideTests(dge01, p.value = 1e-2, lfc=1, p.adjust="BH")
summary(sig)
genes <- rownames(sig)
s.g. <- genes [sig<0 | sig>0]
### decimos que seleccione a los genes GDE entre las 
s.lcpm <- logcounts[s.g., 1:12]
write.csv(s.lcpm,"matriz_heatmap.csv")
Heatmap(as.matrix(s.lcpm))
Heatmap(as.matrix(s.lcpm),row_names_gp = gpar(fontsize=8))
######### PLOT HEATMAP#### CON COLORES
### Recordando que s.lcp es la matriz de log-counts per million (logCPM)
########################################################################
library(ComplexHeatmap)
library(circlize)
### Escalamiento de filas de genes, para que el HeatMap exprese patrones de expresion
### y no de abundancias absolutas!!
### Paso 01
logcounts <- cpm(N01, log=TRUE)
### Selección de genes DEG
### Paso 02
sig <- decideTests(dge01, p.value = 1e-2, lfc = 1, p.adjust.method = "BH")
summary(sig)
## Extacción de genes significativos
genes_sig <- rownames(sig)[which(sig != 0)]
## N° de DEGs o cuantos DEGs hay
length(genes_sig)
### Filtación de logcounts con los genes
s.lcpm <- logcounts[genes_sig, 1:12]
### Escalacion por filas y Visualizar los patrones
mat.scaled <- t(scale(t(s.lcpm)))
## Si se tiene la info de los grupos, creamos anotaciones
### ha_col <- HeatmapAnnotation(Grupo = sample_info$group)
#ha_col <- heatmapAnnotation(Grupo = sample_info$group)
## Heatmap PLOTs
Heatmap(
  mat.scaled,
  name = "logCPM (z-score)",
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 10),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

## Exportamos la matrix, si en el caso es necesario
# write.csv(s.lcpm, "matriz_heatmap.csv")

### PLOT PARA SELECCION DE 100 GENES MAS VARIABLES ##
vars <- apply(s.lcpm, 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:100]
mat.top <- s.lcpm[top_genes, ]
mat.scaled.top <- t(scale(t(mat.top)))

#### HEATMAP PLOT 100 GENES VARIABLES, ANOTACIONES Y TOP GENES
##############################################################
Heatmap(
  mat.scaled.top,
  name = "logCPM",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)
### HEATMAP MAS SEXI!!! MY NIGGAAAAAAAAAAAA!!!! #########
######## GIANMARCO CASTILLO HUACCHO - UNJFSC - TUTORIAL ########

sel_genes <- decideTests(dge01, p.value = 0.01, lfc = 1, p.adjust.method = "BH")
summary(sel_genes)
## 
gene_names <- rownames(sel_genes)
deg_genes <- gene_names[sel_genes != 0]
#### Matriz log-normalizada
logcounts <- cpm(N01, log = TRUE)
deg_matrix <- logcounts[deg_genes, ]
#### Anotaciones
group_labels <- factor(rep(c("GrupoA", "GrupoB"), each = 6))
## DOS GRUPOS
group_colors <- c("GrupoA" = "#1f78b4", "GrupoB" = "#33a02c")

### ESTE ES ANOTACION DEL HEATMAP CON DOS GRUPOS
col_annotation <- HeatmapAnnotation(
  Grupo = group_labels,
  col = list(Grupo = group_colors),
  annotation_name_side = "left"
)
### Heatmap con dos grupos A y B

Heatmap(as.matrix(deg_matrix),
        name = "logCPM",
        top_annotation = col_annotation,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = TRUE)

##### YA TENEMOS EL HEATMAP CON LAS ANOTACIONES
##### AHORA VAMOS HACERLO MAS SEXI, PARA UNA PUBLICACION CIENTIFICA
##### HEATMAP MAS SEXI! MIS CHAM@S!!

h <- Heatmap(
  as.matrix(deg_matrix),
  name = "logCPM",
  col = colorRamp2(
    c(min(deg_matrix), 0, max(deg_matrix)),
    c("#440154", "#FDE725", "#21908C")
  ),
  top_annotation = col_annotation,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = TRUE
)

h

ggsave(plot = h, "Heatmap_Sexi_FINISH.pdf", dpi = 500, width = 18.5, height = 12.9)


#################
######## Introducción al volcano plot ######

head(tags01$table)
plot(tags01$table$logFC, -10*log (tags01$table$PValue))
table<- tags01$table

table$sig <- ifelse(table$FDR<0.01 & table$logFC<0 , "down", "no")
summary(sig)

ts <- data.frame(geneID = rownames(sig), SIG = sig)
head(ts)
names(ts) <- c("geneID","SIG")
head(ts)
table$geneID <- rownames(table)
head(table)

table <- merge(table, ts, by="geneID")
head(table)
library(ggplot2)
library(ggrepel)
ggplot(table, aes(x=logFC, y= -10*log(PValue), col=as.factor(SIG))) + geom_point() +
  theme_classic()

#####################################################################
#####################################################################
#####################################################################
#####################################################################
## VOLCANO PLOT CON GENES DOWN REGULADOS Y UP REGULADOS CON ETIQUETAS

# Preparar tabla
table <- tags01$table
table$geneID <- rownames(table)

# Clasificar genes en up, down y no significativos
table$SIG <- ifelse(table$FDR < 0.01 & table$logFC > 0, "up",
                    ifelse(table$FDR < 0.01 & table$logFC < 0, "down", "no"))

# Seleccionar genes más significativos para anotar (top 10 up y top 10 down)
top_up <- table[table$SIG == "up", ]
top_up <- top_up[order(top_up$FDR), ][1:10, ]

top_down <- table[table$SIG == "down", ]
top_down <- top_down[order(top_down$FDR), ][1:10, ]

genes_to_label <- rbind(top_up, top_down)

# Volcano plot con colores personalizados y anotaciones
ggplot(table, aes(x = logFC, y = -10 * log10(PValue), col = SIG)) +
  geom_point(alpha = 0.7, size = 2) +
  # Colores personalizados
  scale_color_manual(values = c("up" = "#E76F51", "down" = "#457B9D", "no" = "#BDBDBD")) +
  # Anotaciones
  geom_text_repel(data = genes_to_label, aes(label = geneID), size = 3, max.overlaps = 20) +
  # Líneas de referencia
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -10 * log10(0.01), linetype = "dotted", color = "black") +
  # Tema y títulos
  theme_classic(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significance"
  )


###############################################################
###############################################################
#### VOLCANO-PLOT CON COLORES - ETIQUETAS DOWN GENES Y UP GENES

# Obtener la tabla
table <- tags01$table
table$geneID <- rownames(table)

# Clasificar genes como up / down / no significativos
table$SIG <- ifelse(table$FDR < 0.01 & table$logFC > 1, "up",
                    ifelse(table$FDR < 0.01 & table$logFC < -1, "down", "no"))

# Seleccionar genes a anotar (top 10 up y top 10 down por FDR)
top_up <- table[table$SIG == "up", ]
top_up <- top_up[order(top_up$FDR), ][1:10, ]

top_down <- table[table$SIG == "down", ]
top_down <- top_down[order(top_down$FDR), ][1:10, ]

genes_to_label <- rbind(top_up, top_down)


##### Volcano plot con mejoras
ggplot(table, aes(x = logFC, y = -10 * log10(PValue), col =as.factor(SIG))) +
  geom_point(alpha = 0.8, size = 2.5) +
  # Colores personalizados (distintos a rojo/azul)
  scale_color_manual(values = c(
    "up" = "#E76F51",
    "down" = "#3A86FF",
    "no" = "#BDBDBD"
  )) +
  # Anotar genes con ggrepel
  geom_text_repel(data = genes_to_label, aes(label = geneID), size = 3, max.overlaps = 20) +
  # Líneas de referencia para logFC y p-valor
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "purple") +
  geom_hline(yintercept = -10 * log10(0.01), linetype = "dashed", color = "purple") +
  # Temas y etiquetas
  theme_classic(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Regulación"
  )

### VOLCANO PLOT MAS SEXI!########
### PARA PUBLICACIÓN CIENTIFICA###
### UNJFSC - HUACHO CITY!!!#######

ggplot(table, aes(x = logFC, y = -10 * log10(PValue), col = as.factor(SIG))) +
  geom_point(alpha = 0.6, size = 1.8) +  # puntos más pequeños y semitransparentes
  scale_color_manual(values = c(
    "up" = "#E76F51",
    "down" = "#3A86FF",
    "no" = "#BDBDBD"
  )) +
  geom_text_repel(
    data = genes_to_label,
    aes(label = geneID),
    size = 2.8,
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = 15,
    segment.color = "grey70",
    segment.size = 0.2
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -10 * log10(0.01), linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_minimal(base_size = 11) +  # tema más limpio y moderno
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(
    title = "Volcano Plot",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10](p~value)),
    color = "Regulación"
  )

ggsave("Volcano_SEXI_GAAAA.pdf", dpi = 500, width = 8.5, height = 6.9)


############################
############################ ONTOLOGIA GENICA
############################

#install.packages ("installr")
#updateR()         ESTO ES EN WINDOWS PARA LINUX UTILIZAR SUDO

#BiocManager::install("topGO")
library(topGO)
setwd("/home/marcos/Descargas/RNA-seq/RNAseq/")
#leer el archivo de anotaci?n funcional
geneID2GO<-readMappings(file="Cquinoa.GO.tsv")

#ir al archivo Cquinoa.GO, buscar -RA y reemplazar por nada, para eliminarlo
geneID2GO<-readMappings(file="Cquinoa.GO.tsv")
str(head(geneID2GO))
#### visualizacion de los Genes ontologos
View(geneID2GO)
#### Nombre de los genes ontologos
names(geneID2GO)
#archivo de expresi?n diferencial
head("tags01$table")
#seledcciono la tabla y le digo que sea menor o igual a 1 y el logFC 1, 
#por lo menos 1 tratamiento sea el doble del otro
#rownames le pide solo los nombres de los genes, en este ejemplo son 75 genes
rownames(tags01$table[tags01$table$FDR<=0.01 & abs(tags01$table$logFC)>=1, ])


#corro con tags01$table
selectgenes<-rownames (tags01$table[tags01$table$FDR<=0.01 & abs(tags01$table$logFC)>=1, ])
## Nombre de los genes 
geneNames<-names(geneID2GO)
## Visualizacion de los genes
View(geneNames)
 
#no todos los genes seleccionados van a tener una anotaci?n
#geneList son todos los genes que tienen anotaci?n funcional
geneList<- factor(as.integer(geneNames %in% selectgenes))
#genera los nombres
names(geneList) <- geneNames
#geneList solo tiene una lista de 0 y 1 que dice si estaba o no seleccionado
str(geneList)
summary(geneList)
#ontologia de procesos biologicos=BP, el gene List - es la lista de los genes 
#y cuales son seleccionados (1) o 
#no seleccionado (0)
#el GOdata es el que usamos para crear las tablas estadisticas para las sgntes funciones:
### BP: Proceso biologico
### MF: Funcion molecular
#### CC: Componente celular

GOdata<-new("topGOdata", ontology="MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO=geneID2GO)

GOdata
#corremos un test, en el cual no todos los genes tienen funciones
results<-runTest(GOdata, statistic="fisher")
results

#################
#RESULT p < 0.01  es la representaci?n que 43 genes significatios que al azar an tenido una frecuencia 
#la probablilidad de seleccionar aleatoreamente #### quiere decir que la agrupaci?n de datos no fue casualidad

###################
tab<-GenTable(GOdata, classic = results, topNodes=5, numChar=100)
tab
#guardar la tabla
write.table(tab, file="Func.enrich.cq01.DRO.tsv", row.names=F, col.names=T, sep="\t", quote=F)
tab<-GenTable(GOdata, classic=results, topNodes=5)

#tambi?n podemos hacer GODATA para ontologia MF (Funcion molecular)
#qu? tipo de reacci?n bioqu?mica cataliza la prote?na.
GOdatamf<-new("topGOdata", ontology="MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO=geneID2GO)
resultsmf<-runTest(GOdatamf, statistic="fisher")
tab<-GenTable(GOdatamf, classic = resultsmf, topNodes=5, numChar=100)
write.table(tab, file="Func.enrichMF.cq01.DRO.tsv", row.names=F, col.names=T, sep="\t", quote=F)
write.table(tab, file="Func.enrichMF.cq01.DRO.tsv", row.names=F, col.names=T, sep="\t", quote=F)
head(tab)
#si quiero cargar el archivo
#para guardar mi lista de genes
dge <- read.table("CPOSvsDPOS.DGE.tsv", sep="\t", header=T, row.names=1)
rownames(dge[dge$FDR<=0.01 & abs(dge$logFC)>=1, ])



## PARA VISUALIZAR EL ANALISIS DE GO SE NECESITARAN ESTOS PAQUETES
#install.packages("ggplot2")
#install.packages("topGO")
#install.packages("GOstats")
#install.packages("clusterProfiler")

library(GOstats)
library(ggplot2)
library(clusterProfiler)
#BiocManager::install("GOstats")

ID2GO<-readMappings(file="Func.enrichMF.cq01.DRO.tsv")


##VISUALIZACION DE GO CON BARCHART

tab <- GenTable(GOdata, classic = results, topNodes = 10, numChar = 100)
tab$classic <- as.numeric(tab$classic)

tab$Term <- factor(tab$Term, levels = rev(tab$Term))

ggplot(tab, aes(x = Term, y = -log10(classic))) +
  geom_bar(stat = "identity", fill = "#3A86FF") +
  coord_flip() +  # horizontal
  theme_minimal(base_size = 12) +
  labs(
    title = "Enriquecimiento de términos GO (MF)",
    x = "Término GO",
    y = expression(-log[10](p~value))
  )

#### En el resultado anterior se puede observar los barchart del analisis de GO
#### Pero no podemos ver la significancia, a si que aqui mejoramos el plot
#### Por cierto esto esto te sirve para entender de una mejor manera el GO


ggplot(tab, aes(x = Term, y = -log10(classic), fill = Significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "skyblue", high = "darkblue") +
  theme_minimal() +
  labs(
    title = "GO Enrichment (MF)",
    x = "GO term",
    y = expression(-log[10](p~value)),
    fill = "Genes significativos"
  )

ggsave("GO_analisis_Barchart.pdf", dpi = 500, width = 14.5, height = 9.9)


#### GRAFICO DE GO CON DOTPLOT ###
#################################

tab$classic <- as.numeric(tab$classic)
tab$Term <- factor(tab$Term, levels = rev(tab$Term))

ggplot(tab, aes(
  x = Significant,
  y = Term,
  size = Annotated,
  color = -log10(classic)
)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  scale_size(range = c(3, 8)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "GO Enrichment Dotplot (MF)",
    x = "Genes significativos",
    y = "Término GO",
    size = "Genes anotados",
    color = expression(-log[10](p~value))
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("GO_analisis_Dotplot.pdf", dpi = 500, width = 18.5, height = 10.9)



### GRAFICAR GO CON DIAGRAMA DE PUNTOS.
##########################################

tab$Term <- factor(tab$Term, levels = rev(tab$Term))  # Orden de arriba hacia abajo
tab$classic <- as.numeric(tab$classic)  

ggplot(tab, aes(
  x = Significant,
  y = Term,
  size = Annotated,
  color = -log10(classic)
)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "skyblue", high = "darkblue") +
  scale_size(range = c(3, 10)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Gráfico de burbujas - Enriquecimiento GO (MF)",
    x = "Genes significativos",
    y = NULL,
    size = "Genes anotados",
    color = expression(-log[10](p~value))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11)
  )
