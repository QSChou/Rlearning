install.packages("quantmod")
library(quantmod)
hist_draw <- function(obj) {
  barChart(obj["2016-01-01::2016-09-19"])
  chartSeries(obj,subset ="2015-01-20::2016-09-19",theme = "white")
  ## 均線
  #SMA(Cl(obj))
  addTA(SMA(Cl(obj)),on=1,col="blue")
  addTA(SMA(Cl(obj),n=20),on=1,col="red")
  ## MACD
  addMACD()
  ## Bollinger Band
  addBBands()
}
# to get the apple stock
#getSymbols("AAPL")
stk_tw3481 <- get(getSymbols("3481.tw"))
hist_draw(stk_tw3481)
stk_tw4504 <- get(getSymbols("4504.tw"))
hist_draw(stk_tw4504)
##回測前制作業
