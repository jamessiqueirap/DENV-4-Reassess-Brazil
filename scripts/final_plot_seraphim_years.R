library(devtools)
library(seraphim)
library(diagram)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(RColorBrewer)
library(ape)
library(viridis)

#setwd("seraphim/")

#### FUNÇÃO DE EXTRAÇÃO DAS INFORMAÇÕES DA TREE###
mccExtractions = function(mcc_tre, mostRecentSamplingDatum, coordinates_file) {
  # Carregar o CSV de coordenadas
  coordinates_df = read.csv(coordinates_file, stringsAsFactors = FALSE)
  
  # Criar tabela para armazenar dados da árvore
  mcc_tab = matrix(nrow=dim(mcc_tre$edge)[1], ncol=7)
  colnames(mcc_tab) = c("node1", "node2", "length", "startLon", "startLat", "endLon", "endLat")
  mcc_tab[,c("node1","node2")] = mcc_tre$edge
  mcc_tab[,c("length")] = mcc_tre$edge.length
  
  # Loop sobre os nós para preencher as coordenadas de fim (endLon, endLat)
  for (i in 1:length(mcc_tre$annotations)) {
    annotations = mcc_tre$annotations[[i]]
    Region_name = annotations$Region
    
    # Buscar as coordenadas correspondentes no CSV
    coords = coordinates_df[coordinates_df$Region == Region_name, ]
    
    if (nrow(coords) > 0) {
      mcc_tab[i, c("endLon", "endLat")] = c(coords$longitude, coords$latitude)
    } else {
      mcc_tab[i, c("endLon", "endLat")] = NA  # Caso não encontre a divisão, marcar como NA
      print(paste("Divisão não encontrada:", Region_name))  # Imprimir o nome da divisão não encontrada
    }
  }
  
  # Preencher as coordenadas de início (startLon, startLat)
  for (i in 1:length(mcc_tre$annotations)) {
    index = which(mcc_tab[,"node2"] == mcc_tab[i,"node1"])
    if (length(index) > 0) {
      mcc_tab[i,c("startLon","startLat")] = mcc_tab[index,c("endLon","endLat")]
    } else {
      annotations = mcc_tre$root.annotation
      Region_name = annotations$Region
      
      # Buscar as coordenadas para o nó raiz
      coords = coordinates_df[coordinates_df$Region == Region_name, ]
      
      if (nrow(coords) > 0) {
        mcc_tab[i, c("startLon", "startLat")] = c(coords$longitude, coords$latitude)
      } else {
        mcc_tab[i, c("startLon", "startLat")] = NA
        print(paste("Divisão não encontrada para o nó raiz:", Region_name))  # Imprimir divisão não encontrada para o nó raiz
      }
    }
  }
  
  # Restante da função original (cálculo de alturas, anos, etc.)
  l = length(mcc_tab[,1])
  ll = matrix(1:l, nrow=l, ncol=l)
  ll[] = 0
  for (j in 1:l) {
    subMat = mcc_tab[j, 2]
    subMat = subset(mcc_tab, mcc_tab[,2] == subMat)
    ll[j,1] = subMat[,3]
    subMat = subMat[1, 1]
    subMat1 = subset(mcc_tab, mcc_tab[,2] == subMat)
    for (k in 1:l) {
      if (nrow(subMat1) > 0) {
        ll[j,k+1] = subMat1[,3]
        subMat2 = subMat1[1,1]
        subMat1 = subset(mcc_tab, mcc_tab[,2] == subMat2)
      }
    }
  }
  
  endNodeL = rowSums(ll)
  mcc_tab = cbind(mcc_tab, endNodeL)
  startNodeL = matrix(1:l, nrow=l, ncol=1)
  startNodeL[] = 0
  for (j in 1:l) {
    r = mcc_tab[j,1]
    s = subset(mcc_tab, mcc_tab[,2] == r)
    for (k in 1:l) {
      if (nrow(s) > 0) {
        startNodeL[j,1] = s[,8]
      }
    } 
  }
  colnames(startNodeL) = "startNodeL"
  mcc_tab = cbind(mcc_tab, startNodeL)
  
  maxEndLIndice = which.max(mcc_tab[,"endNodeL"])
  maxEndL = mcc_tab[maxEndLIndice,"endNodeL"]
  endYear = matrix(mcc_tab[,"endNodeL"] - maxEndL)
  endYear = matrix(mostRecentSamplingDatum + (endYear[,1]))
  startYear = matrix(mcc_tab[,"startNodeL"] - maxEndL)
  startYear = matrix(mostRecentSamplingDatum + (startYear[,1]))
  colnames(startYear) = "startYear"
  colnames(endYear) = "endYear"
  
  mcc_tab = cbind(mcc_tab, startYear, endYear)
  mcc_tab = mcc_tab[order(mcc_tab[,"startYear"], decreasing=F),]
  
  mcc_tab1 = mcc_tab[1,]
  mcc_tab2 = mcc_tab[2:dim(mcc_tab)[1],]
  mcc_tab2 = mcc_tab2[order(mcc_tab2[,"endYear"], decreasing=F),]
  mcc_tab = rbind(mcc_tab1, mcc_tab2)
  
  return(mcc_tab)
}

# 1. Extraindo informações da árvore MCC

mcc_tre = readAnnotatedNexus("Run_07_SkyGrid_MCMC_04_50mi")
mostRecentSamplingDatum = 2024.7814
mcc_tab = mccExtractions(mcc_tre, mostRecentSamplingDatum, "Region_trait_Coordenates.csv")
write.csv(mcc_tab, "denv_cood.csv", row.names=FALSE, quote=FALSE)
coordinates_df = read.csv("Region_trait_Coordenates.csv", stringsAsFactors = FALSE)

# Determinar os índices correspondentes ao primeiro e último ano de transmissão
minYear <- min(mcc_tab[,"endYear"], na.rm = TRUE)
maxYear <- max(mcc_tab[,"endYear"], na.rm = TRUE)

# Ajustar o intervalo para começar em 1976 e terminar em 2020
first_year <- max(1975, floor(minYear))  # Garante 1976 como limite inferior
last_year <- min(2020, ceiling(maxYear)) # Limita o último ano a 2020

# Criar uma paleta de cores contínua para o intervalo de anos ajustado
colour_scale <- colorRampPalette(brewer.pal(11, "RdYlBu"))(last_year - first_year + 1)

# Ajustar os índices para o intervalo ajustado
endYears_indices <- pmax(1, pmin(round((mcc_tab[,"endYear"] - first_year) + 1), length(colour_scale)))

# Mapear as cores usando os índices ajustados
endYears_colours <- colour_scale[endYears_indices]

# 3. Plotando a árvore MCC
pdf("new_alta_qualidade_delimitado.pdf", width = 16, height = 8, onefile = TRUE, useDingbats = FALSE)

# Carregar fronteiras nacionais
borders <- ne_countries(scale = "large", returnclass = "sf")
borders <- st_transform(borders, crs = 4326) # Certificando-se de usar o CRS correto

# Carregar fronteiras estaduais (subdivisões) e filtrar para o Brasil
states <- ne_states(country = "Brazil", returnclass = "sf")
states <- st_transform(states, crs = 4326) # Certificando-se de usar o CRS correto

# Criar um objeto sf e obter o bbox
sf_data <- st_as_sf(data.frame(lon = c(-180, 180), lat = c(-90, 90)), coords = c("lon", "lat"), crs = 4326)
bbox <- st_bbox(sf_data)

# Plotando usando base graphics
par(mar = c(0, 0, 0, 0), oma = c(1.2, 3.5, 1, 0), mgp = c(0, 0.4, 0), lwd = 0.2, bty = "o")
plot(st_geometry(borders), col = "white", border = "gray50", lwd = 0.5)

# Adicionando as fronteiras estaduais do Brasil
plot(st_geometry(states), col = "white", border = "gray50", lwd = 0.5, add = TRUE)

valid_coords <- function(lon, lat) {
  all(lon >= -180 & lon <= 180 & lat >= -90 & lat <= 90, na.rm = TRUE)
}

for (i in 1:nrow(mcc_tab)) {
  if (valid_coords(mcc_tab[i, "startLon"], mcc_tab[i, "startLat"]) &&
      valid_coords(mcc_tab[i, "endLon"], mcc_tab[i, "endLat"])) {
    curvedarrow(
      cbind(mcc_tab[i, "startLon"], mcc_tab[i, "startLat"]),
      cbind(mcc_tab[i, "endLon"], mcc_tab[i, "endLat"]),
      arr.length = 0, arr.width = 0, lwd = 2, lty = 1, 
      lcol = endYears_colours[i], arr.col = NA, arr.pos = FALSE, curve = 0.3
    )
  }
}

# Mesclando os dados de mcc_tab com o arquivo de coordenadas
merged_data <- merge(
  mcc_tab,
  coordinates_df,
  by.x = "startLon",
  by.y = "longitude",  # Substitua 'longitude' pelo nome correto da coluna em coordinates_df
  all.x = TRUE
)

# Plotando pontos e nomes das localidades
for (i in seq_len(nrow(merged_data))) {
  if (valid_coords(merged_data[i, "endLon"], merged_data[i, "endLat"])) {
    points(merged_data[i, "endLon"], merged_data[i, "endLat"], pch = 16, col = endYears_colours[i], cex = 1)
    points(merged_data[i, "endLon"], merged_data[i, "endLat"], pch = 1, col = "gray10", cex = 1)
  }
}

# Usar os valores do bbox para plotar o retângulo
rect(
  xleft = bbox["xmin"],
  ybottom = bbox["ymin"],
  xright = bbox["xmax"],
  ytop = bbox["ymax"],
  xpd = TRUE, lwd = 0.2
)

# Gerar sequência de anos inteiros para a escala contínua
years_seq <- seq(first_year, last_year, by = 1)

# Criar um raster representando os anos (escala contínua)
rast <- raster(matrix(years_seq, nrow = 1, ncol = length(years_seq)))
extent(rast) <- c(first_year, last_year, 0, 1)

# Plotar a escala de cores ajustada (gradiente contínuo)
plot(
  rast, legend.only = TRUE, add = TRUE, col = colour_scale, 
  legend.width = 1, legend.shrink = 0.75, horizontal = TRUE,
  smallplot = c(0.2, 0.80, 0.02, 0.03), # Ajuste da posição da legenda
  axis.args = list(
    at = seq(first_year, last_year, by = 5), # Ticks a cada 5 anos
    labels = seq(first_year, last_year, by = 5), # Rótulos dos anos
    cex.axis = 1, tck = -0.7, col.axis = "black",
    line = -0.020, lwd = 0.5
  ),
  legend.args = list(
    text = "Years", side = 3, cex = 1, line = 0.5, col = "black"
  )
)

dev.off()

