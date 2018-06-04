library(coinmarketcapr)
library(treemap)
library(ggplot2)
library(FactoMineR)
library(stats)
library(factoextra)
#plot_top_5_currencies()

market_today <- get_marketcap_ticker_all()
head(market_today[,1:8])

df1 <- na.omit(market_today[,c('symbol','market_cap_usd', 'percent_change_24h', 'price_btc', 'percent_change_7d', 'total_supply','X24h_volume_usd', 'available_supply')])

#as numeric
df1$market_cap_usd <- as.numeric(df1$market_cap_usd)
df1$percent_change_24h <- as.numeric(df1$percent_change_24h)
df1$price_btc <- as.numeric(df1$price_btc)
df1$percent_change_7d <- as.numeric(df1$percent_change_7d)
df1$total_supply <- as.numeric(df1$total_supply)
df1$X24h_volume_usd <- as.numeric(df1$X24h_volume_usd)
df1$available_supply <- as.numeric(df1$available_supply)

#df1$formatted_market_cap <-  paste0(df1$id,'\n','$',format(df1$market_cap_usd,big.mark = ',',scientific = F, trim = T))

cryptolist <- df1[1:50, 2:8]
row.names(cryptolist) <- df1[1:50,1]

#Tree Map
#treemap(df1, index = 'formatted_market_cap', vSize = 'market_cap_usd', title = 'Cryptocurrency Market Cap', fontsize.labels=c(12, 8), palette='RdYlGn')

res.pca <- prcomp(cryptolist, scale = TRUE)
fviz_eig(res.pca)

#graph of coins, similar group together
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             axes = c(1,3)
)
#graph variables - positively correlated point to same side. negative point opposite
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE ,    # Avoid text overlapping
             axes = c(1,3)
)

#Biplot
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  ,# Individuals color
                axes = c(1,2)
)
