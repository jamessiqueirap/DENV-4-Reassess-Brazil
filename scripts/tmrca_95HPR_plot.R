# Carregar o ggplot2 e dplyr
library(ggplot2)
library(dplyr)

# Função para converter ano decimal para data
decimal_year_to_date <- function(year_decimal) {
  # Extrair o ano
  year <- floor(year_decimal)
  # Extrair a fração do ano
  fraction <- year_decimal - year
  # Calcular o número de dias no ano (para anos bissextos e comuns)
  days_in_year <- ifelse((year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0), 366, 365)
  # Converter fração para dias
  date <- as.Date(paste0(year, "-01-01")) + round(fraction * days_in_year)
  return(date)
}

# Ler os dados do arquivo TSV
df <- read.csv("hpd_95%.tsv", sep="\t", header=TRUE)

# Converter apenas as colunas numéricas (tmrca, hpd_low, hpd_high)
# Identificar as colunas específicas por nome ou por índice
df$tmrca_value <- as.numeric(gsub(",", ".", df$tmrca_value))
df$low_value <- as.numeric(gsub(",", ".", df$low_value))
df$high_value <- as.numeric(gsub(",", ".", df$high_value))

# Aplicar a conversão de ano decimal para data
df$tmrca_value <- decimal_year_to_date(df$tmrca_value)
df$low_value <- decimal_year_to_date(df$low_value)
df$high_value <- decimal_year_to_date(df$high_value)

# Ajustar a localidade para garantir rótulos em inglês
Sys.setlocale("LC_TIME", "C")

# Ordenar as divisões de modo decrescente por tmrca_value
df <- df %>%
  arrange(desc(tmrca_value)) %>%
  mutate(division = factor(division, levels = unique(division)))

# Criar o gráfico
p <- ggplot(df, aes(y = division)) + 
  geom_errorbarh(aes(xmin = low_value, xmax = high_value), height = 0.2, color = "black") +
  geom_segment(aes(x = tmrca_value - 30, xend = tmrca_value + 30, y = division, yend = division), 
               size = 0.5, color = "white") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    text = element_text(size = 12, color = "black"), # Texto geral
    axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5), # Alinhamento dos rótulos do eixo X
    axis.text.y = element_text(color = "black"), # Texto do eixo Y normal
    axis.title = element_text(color = "black"), # Títulos dos eixos
    plot.title = element_text(color = "black"), # Título do gráfico
    panel.grid.minor = element_blank(), # Remove grades menores
    panel.grid.major.y = element_blank(), # Remove grades principais horizontais
    panel.grid.major.x = element_blank(), # Remove grades principais verticais
    axis.line = element_line(color = "black"), # Adiciona linhas nos eixos X e Y
    axis.ticks = element_line(color = "black", size = 0.5) # Adiciona tracinhos nos eixos
  )

print(p)
# Salvar o gráfico em PDF
ggsave("tmrca_and_952.pdf", plot = p, width = 14, height = 7, dpi = 400)
