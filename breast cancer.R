library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(caret)
library(psych)
library(pastecs)
library(mice)
library("corrplot")
library(factoextra)
library(ggrepel)
library(clustertend)
library(fpc)
library(tidyverse)
library(cluster)
#### VERİ SETİ BİLGİSİ ####
getwd()
setwd("C:/Users/gamze/OneDrive/Masaüstü")
list.files()
df<-read.csv("wdbc.data"      ,header=TRUE)
head(df)
View(df)
colnames(df)=c("ID","diagnosis","radius","texture","perimeter",
"area","smoothness","compactness","concavity","concave_points"
,"symmetry","fractal_dimension")
df<-df[,-c(1,2,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)]
df
head(df)
str(df)

## Tanımlayıcı İstatistikler


summary(df)
head(df)
stat.desc(df)

apply(df,2,sd)
apply(df,2,var)

### **Kutu Grafikleri**

boxplot(df, Vertical=T)

boxplot(df$radius, horizantal = T, xlab= "radius")
boxplot(df$ texture, horizantal = T,xlab =" texture")
boxplot(df$perimeter , horizantal = T,xlab= "perimeter ")
boxplot(df$ area , horizantal = T,xlab= " area ")
boxplot(df$smoothness, horizantal = T, xlab="smoothness")
boxplot(df$compactness, horizantal = T,xlab = "compactness")
boxplot(df$concavity, horizantal = T, xlab="concavity")
boxplot(df$concave_points , horizantal = T,xlab="concave_points ")
boxplot(df$symmetry, horizantal = T,xlab="symmetry")
boxplot(df$fractal_dimension, horizantal = T, xlab="fractal_dimension")

attach(df)

### **Değişkenler Arasındaki Korelasyon**

corr <- cor((df), method = "pearson")
corr

corrplot.mixed(corr, lower="pie",upper="number")  


### ***Hopkins Testi***
set.seed(123)
hopkins.data <- hopkins(df,n=nrow(df)-1)
hopkins.data

dfScaled <- preProcess(df, method=c("center", "scale"))
dfScaled <- predict(dfScaled, newdata = df)
boxplot(dfScaled)

## Temel Bileşenler Analizi

##korelasyon matrisi üzerinden
dfScaled.pca <- prcomp(dfScaled, center = TRUE, scale. = TRUE)  
summary(dfScaled.pca)
### varyans kovaryans matrsisi üzerinden
dfscaled.pca<-prcomp(dfScaled)
summary(dfScaled.pca)

(dfScaled.pca$sdev)^2

fviz_eig(dfScaled.pca)

x <- fa.parallel(dfScaled, fm="pa", fa="both", n.iter=1)

##ÖZVEKTÖRLER##
dfScaled.pca$rotation[,1:2] 
###yeni gözlem değerleri
dfscaled.pca$x
#iki bileşene ait yeni veri seti
df_new <- dfScaled.pca$x[,1:2]
df_new
##DEĞİŞKENLERİN KATKILARI###
fviz_contrib(dfScaled.pca, choice = "var", axes = 1,top = 10)
fviz_contrib(dfScaled.pca, choice = "var", axes = 2,top = 10)
####PC1-PC2
fviz_pca_ind(dfScaled.pca, axes = c(1,2),
             col.ind = "cos2",gradient.cols = c("#00AFBB","#E7B800","#FC4E07"),
             repel = TRUE )

df[c(78,461,504,212),]
summary(df)

fviz_pca_var(dfScaled.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)


res.var <- get_pca_var(dfscaled.pca)
res.var$cos2 
res.var$contrib 
corrplot(res.var$cos2, is.corr=FALSE) ##bunu ekle çok fazla değişken var
corr
cor(dfScaled.pca$x[,1],dfScaled.pca$x[,2])
# Results for individuals
res.ind <- get_pca_ind(dfscaled.pca)
res.ind$coord [,1:2]          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation

###UZAKLIK ÖLÇÜLERİ###

dist.euc=get_dist(df_new, method="euclidean")
fviz_dist(dist.euc)

##MOR RNKLER UZAKLIĞIN EN FAZLA OLDUĞU YERLER



dist.cor=get_dist(df, method="pearson")
fviz_dist(dist.cor)

###KÜMELEME ANALİZİ###
## K MEANS##
###ELBOW##
set.seed(123)
fviz_nbclust(df_new,kmeans,method = "wss",nstart = 25)
##SILIHUTTE##
set.seed(123)
fviz_nbclust(df_new,kmeans,method = "silhouette") #for average silhouette width
##GAP İSTATİSTİĞİ##
set.seed(123)
fviz_nbclust(df_new,kmeans,method = "gap_stat", nboot=50)

##K MEANS##
set.seed(123)
km_res <- kmeans(df_new, 2, nstart=25) 
print(km_res)
fviz_cluster(km_res, data = df_new,
             ellipse.type = "convex",         # Concentration ellipse
             star.plot = TRUE,            # Add segments from centroids to items
             repel = TRUE,               # Avoid label overplotting (slow)
             ggtheme = theme_minimal())




set.seed(123)
km_res <- kmeans(df_new, 6, nstart=25) 
#print(km_res)


fviz_cluster(km_res, data = df_new,axes = c(1,2),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)


set.seed(123)
km_res <- kmeans(df_new, 3, nstart=25) 
#print(km_res)


fviz_cluster(km_res, data = df_new,axes = c(1,2),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)

##K MEDOIDS##
#  **Küme Sayısının Belirlenmesi**

set.seed(123)
fviz_nbclust(df_new,pam,method = "wss") ##4 
fviz_nbclust(df_new,pam, method = "silhouette") ##2VE5##

set.seed(123)
pam_df_2 <- pam(df_new,4)
print(pam_df_2)
fviz_cluster(pam_df_2,
             ellipse.type = "convex"
             , repel = TRUE
             , axes = c(1,2))
#k-medoids pca
set.seed(123)
pam_df_2 <- pam(df_new,4)
fviz_cluster(pam_df_2,
             ellipse.type = "convex"
             , repel = TRUE
             , axes = c(1,2)

set.seed(123)
pam_df_2 <- pam(df_new,2)
print(pam_df_2)
fviz_cluster(pam_df_2,
             ellipse.type = "convex"
             , repel = TRUE
             , axes = c(1,2)
             
)

set.seed(123)
pam_df_2 <- pam(df_new,5)
print(pam_df_2)
fviz_cluster(pam_df_2,
             ellipse.type = "convex"
             , repel = TRUE
             , axes = c(1,2)
             
)

###CLARA###
#**Küme Sayısının Belirlenmesi**

set.seed(123)
fviz_nbclust(df_new,clara,method = "wss") ## 3 ya ada 5
fviz_nbclust(df_new,clara, method = "silhouette") ## 2 ve 4

set.seed(123)
print(df_clara)

df_clara <- clara(df_new,3, samples = 50, pamLike = TRUE)

fviz_cluster(df_clara,
             ellipse.type = "convex",geom = "point"
             , pointsize = 1,
             ggtheme = theme_classic()
             
             ,axes = c(1,2))                             ###en uygunu bu###


set.seed(123)
df_clara <- clara(df_new,5, samples = 50, pamLike = TRUE)
#print(df_clara)

fviz_cluster(df_clara,
             ellipse.type = "t",geom = "point"
             , pointsize = 1,
             ggtheme = theme_classic()
             
             ,axes = c(1,2)     
)

set.seed(123)
df_clara <- clara(df_new,2, samples = 50, pamLike = TRUE)
print(df_clara)

fviz_cluster(df_clara,
             ellipse.type = "t",geom = "point"
             , pointsize = 1,
             ggtheme = theme_classic()
             
             ,axes = c(1,2)     
)

set.seed(123)
df_clara <- clara(df_new,4, samples = 50, pamLike = TRUE)
print(df_clara)

fviz_cluster(df_clara,
             ellipse.type = "t",geom = "point"
             , pointsize = 1,
             ggtheme = theme_classic()
             
             ,axes = c(1,2)     
)
###SONUÇ##
set.seed(123)
df_clara <- clara(df_new,3, samples = 50, pamLike = TRUE)
print(df_clara)


set.seed(123)
df_clara <- clara(df_new, 3,  samples = 50, pamLike = TRUE)
print(df_clara)





fviz_cluster(df_clara, data = df_new,
             ellipse.type = "convex", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)
##kümeleme sonuçları##
aggregate(df, by=list(cluster=df_clara$cluster), mean)

summary(df)