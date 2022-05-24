########################################	Proyecto enfermedades monogénicas	################################################
#																															   #
#Autor: Luis Enrique Ramírez Serrano                                                                                           #
#																						   									   #
#Fecha de última actualización: 12/05/2022                                                                                     #
#																					                                           #
#Resumen de funcionalidad: En este programa se termina de limpiar la base de datos de clinvar, pues hay algunas variantes ya   #
#no registradas en la versión del genoma GRCh38. Se lee y limpia la anotación del genoma correspondiente a la versión GRCh38.  #
#Se transforma la base de datos de clinvar y la anotación del genoma de referencia en objetos del paquete GRanges para poder   #
#mapear las variantes y posteriormente se tabulan y grafican las proporciones en los distintos biotipos regristrados.          #
#                                                  																			   #
#Dependencias:                                                                                                                 #
#	clinvar_simp.txt: Forma simplificada y curada de el resultado de la búsqueda de unión entre OMIM y clinvar, proveniente del#
#					  archivo clinvar result.txt                                                                               #
#	GRCh38_latest_genomic.gff: última versión de la anotación del genoma humano, correspondiente a la versión GRCh38           #
#Output: gráfica con los porcentajes de los biotipos y tabla con la concurrencia                                               #
################################################################################################################################

##cut -f1,2,10,11 clinvar_result.txt | grep -E .+$'\t'.+$'\t'.+'\|' -v > clinvar_simp.txt## 
##Esta linea sirve para procesar el archivo clinvar_result.txt, pues solo necesitamos las pociciones genómicas de la versión GRch38, además tiene redundancias en el nombre del cromosoma##
##También se eliminaron manualmente algunas líneas que no contienen información cromosómica

#Se llama a las librerías correspondientes##
library(Biobase)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(GenomicAlignments)
library(ggpercent_plot)
library(stringr)
#El archivo clinvar_simp.txt(puede editarse la ubicación del archivo) es leido y se le quitan las línes que no tengan un valor en locación, esto porque no están registradas esas variantes en la varsión GRCh38 del genoma de referencia
clinvar<-read.table('./clinvar_simp.txt',header=T, sep = "\t")
clinvar<-clinvar[!is.na(as.numeric(clinvar$GRCh38Location)),]

#Lista que contiene los nombres de todos los cromosomas en formato de NCBI RefSeq
ChromosomeNCBI<-c('NC_000001.11','NC_000002.12','NC_000003.12','NC_000004.12','NC_000005.10','NC_000006.12','NC_000007.14','NC_000008.11','NC_000009.12','NC_000010.11','NC_000011.10','NC_000012.12','NC_000013.11','NC_000014.9','NC_000015.10','NC_000016.10','NC_000017.11','NC_000018.10','NC_000019.10','NC_000020.11','NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10','NC_012920.1')
#Como nombres se ponen los números de los cromosomas para poder sustituirlos
names(ChromosomeNCBI)<-unique(clinvar$GRCh38Chromosome)

##Reemplaza los números de los cromosomas por su nombre en NCBI RefSeq
clinvar$GRCh38Chromosome<-ChromosomeNCBI[clinvar$GRCh38Chromosome]

#Creación un objeto de rangos genómicos a partir de las variaciones encontradas en OMIM y clinvar
mendeliandsGR<-GRanges(seqnames = clinvar$GRCh38Chromosome,ranges = IRanges(start = as.numeric(clinvar$GRCh38Location),end=as.numeric(clinvar$GRCh38Location)))

#Lectura de la anotación del genoma de referencia (puede editarse la dirección del archivo)
reference<-read.table('./GRCh38_latest_genomic.gff', sep="\t")
#Solo se utilizan las primeras 7 líneas correspondientes al cromosoma, [], biotipo, inicio, fin, . y cadena
reference<-unique(reference[,1:7])
#Los biotipos correspondientes a "region" y "match" son redundantes, pues representan a los cromosomas o regiones genómicas y estos biotipos no nos interesan
reference<-reference[!(reference$V3 %in% c("region", "match")),]
#Para simplificar la transformación de la anotación del genoma de referencia se nombran cada una de sus columnas
colnames(reference)<-c("Seqnames", "Source", "type", "Start", "End", ".", "Strand")
#creación del objeto Granges con el genoma de referencia
genome_ref<-as(reference, "GRanges")

#Conteo de los sobrelapes para saber cuantas variantes mapearon y cuantas no
sum(overlapsAny(mendeliandsGR, genome_ref))
sum(!overlapsAny(mendeliandsGR, genome_ref))

#Con findOverlaps Obtenemos los índices de los sobrelapes en una matriz
overlaps<-as.matrix(findOverlaps(mendeliandsGR, genome_ref))


#Con la matríz de sobrelapes que tiene los índices podemos buscarlos en el genoma de referencia y colocarlos en la lista de biotipos a los que mapea cada variante (ya que pueden ser varios se crea una lista)
temp<-unique(overlaps[,1])
biotype<-vector(mode="list", length=length(temp))
for(i in temp){
	biotype[[i]]<-as.vector(unique(genome_ref[overlaps[overlaps[,1]==i,2],]$type))
}
#Lista temporal para guardar los biotipos ya que la lista será modificada mas adelante
biotypos2<-biotype

#Debido a que no existe el biotipo de "intrón" consideramos a todo lo que contenga el biotipo "gene", "mRNA" y "transcript" como un intrón si no contiene alguna otra notación (Como CDC o exón)
categories<-c("gene","mRNA","transcript")
for(i in temp){ 
	if (length(biotypos2[[i]])==1) {
		if(biotypos2[[i]]=="CDS"){ 
			biotypos2[i]<-"exon"
		}
	}
	if (all(biotypos2[[i]] %in% categories)) { 
		biotypos2[[i]]<-c(biotypos2[[i]],"intron")
	}
}

#No todas las notaciones contienen biotipo, por lo cual pueden quedar algunos Null que son eliminados
biotypos2<-biotypos2[!(sapply(biotypos2, is.null))]

#Al ya haber definido lo que es un intrón y un exón los biotipos "gene", "mRNA", "transcript" y "CDS" se vuelven redundantes, por lo que es necesario excluirlos del análisis.
redundant<-c("gene","mRNA","transcript","CDS")
table<-table(unlist(biotypos2))
table<-table[!(names(table) %in% redundant)]

#Guardado de la tabla con las frecuencias y concurrencias de cada biotipo
biological<-names(table)
conc<-unname(table)
percent<-conc/sum(conc)*100
table_final<-data.frame(biological, as.vector(conc), as.vector(percent), stringsAsFactors=F)
colnames(table_final)<-c("biotipo", "concurrencia", "porcentaje")
table_final<-table_final[rev(order(table_final$porcentaje)),]

write.table(table_final, "Biotipos.txt", row.names=F, quote=F, sep="\t")

#Gráfica de los porcentajes, si el porcentaje es menor al 0.2% se unirán en la categoria "other, para no imprimir todas las opciones"
data_plot<-(NULL)
data_plot$percent<-c(table_final$porcentaje[table_final$porcentaje>0.2], sum(table_final$porcentaje[table_final$porcentaje<0.2]))
data_plot$names<-c(table_final$biotipo[1:length(data_plot$percent)-1],"other")
data_plot<-as.data.frame(data_plot)
data_plot$names<-paste(data_plot$names, paste(round(data_plot$percent,2), "%", sep=""), sep=": ")
data_plot$names<-factor(data_plot$names, levels=data_plot$names)

colors<-c("#E36B2C", "#BF80FF", "#CC3300", "#8C4966", "#00C0B8", "#00BF7D","#72B000", "#BB9D00", "#9933FF","#024a86")
percent_plot<-ggplot(data_plot, aes(x="", y=percent, fill=(names))) +
	geom_bar(width=1, stat="identity", position="stack")+
	coord_polar("y") +
	labs(x=NULL, y=NULL, fill=NULL, title="Porcentaje de mutaciones asociadas a
		distintos biotipos") +
	theme(axis.line=element_blank(),
		axis.text=element_blank(),
		axis.ticks=element_blank()) +
	scale_fill_manual(values=colors) +
	theme(legend.text=element_text(size=9),
		plot.title=element_text(size=17))
percent_plot
ggsave("Gráfica1.jpg", percent_plot, width=7, height=8.4, units="in", scale=1)
