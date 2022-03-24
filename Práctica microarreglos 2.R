### Práctica expresión diferencial y anotación
library(Biobase)
library(limma)

#Cargamos la tabla
exprs = read.table("expr_normalizada.txt", header=TRUE, row.names=1)

#Revisar el contenido
head(exprs)

#Boxplot de la expresión de genes diferenciados con color

boxplot(exprs, col=rainbow(6))

#Vamos a diferenciar las muestras KO y las WT

types = factor(c("KO","KO","KO","WT","WT","WT"))
types

### hacemos una matriz

design = model.matrix(~ 0+types)
colnames(design) = levels(types)
design

#hacemos una comparación de las muestras 

contMatrix = makeContrasts(KO-WT, levels=design)
contMatrix
help("makeContrasts")


#ajustamos nuestros datos a un modelo lineal

fit = lmFit(exprs,design)

#estimamos los contrastes 

fit2 = contrasts.fit(fit,contMatrix)

#comprimimos las varianzasa un valor lineal

fit2 = eBayes(fit2)

#tomamos los  con genes (sondas) con mejor valor de p

topTable(fit2, number=20, sort.by="p")

#### descubrimos la identidad de los genes

library(mouse4302.db)
library(AnnotationDbi)
library(stats4)
library(IRanges)
library(S4Vectors)

#observamos los datos que contiene
mouse4302()

#fit2 contiene la lista de genes que nos interesa
#extraemos los probes de fit 2

probes = fit2$genes$ID

#Con mouse extraemos los nombres descriptivos

description = mget (probes, mouse4302GENENAME)

# Nombre común

symbols = mget(probes,mouse4302SYMBOL)

#Identificador común

entrezids = mget(probes,mouse4302ENTREZID)

#ponemos estos datos en fit2
fit2$genes$EntrezID = unlist(entrezids)

fit2$genes$Symbol = unlist(symbols)

fit2$genes$Description = unlist(description)

head(fit2$genes)

#volvemos a visualizar los genes con valores de p más significativos

topTable(fit2, number=20, sort.by="p")

#podemos generar tablas con los valores que deseemos
#fold chznge mayor a 1.5
#p value menor a 0.5

deTable = topTable(fit2, number=nrow(fit2), lfc=log2(1.5), p.value=0.05)
dim(deTable)
deTable

#tabla con los datos del microarreglo ordenada por el logFC

fullTable = topTable(fit2, number=nrow(fit2), sort.by="logFC")
dim(fullTable)
### Graficamos un volcano plot para identificar mejor los genes de importancia

volcanoplot(fit2, highlight=10, names=fit2$genes$Symbol)

# ¿Cuántas sondas predices como diferencialmente expresadas?
# 219

# ¿Cuántas decrecen y cuántas aumentan su expresión en el KO?

# ¿Cuántos genes únicos hay en estas listas
unicos <- unique(fit2)
length(unicos)
#24

#Genera una figura de volcán manualmente, comparando explicitamente logFC contra -log10(P.Value) e incluyendo todas las sondas

volcanoplot(fit2, xlab = "Log2 Fold Change", ylab = "p-value", highlight=10, names=fit2$genes$Symbol)

#Colorea de rojo las sondas que encontramos como diferencialmente expresadas.

volcanoplot(fit2, xlab = "Log2 Fold Change", ylab = "p-value", highlight=10, names=fit2$genes$Symbol, hl.col=red)

