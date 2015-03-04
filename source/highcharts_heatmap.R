library(rCharts)
library(rjson)
library(dplyr)

#modified from https://github.com/stefanwilhelm/heatmap. thx to Stefan Wilhelm 

highcharts.heatmap <- function(plot.data, dataType, showLabels=F, useWithShiny=T){
  
  # create Highcharts object
  p <- Highcharts$new()
  
  # use type='heatmap' for heat maps
  p$chart(zoomType = "xy", type = 'heatmap')
  #p$title(text='Sales per employee per weekday')
  p$xAxis(min = 1, max = max(plot.data$x, na.rm=T))
  p$yAxis(min = 1, max = max(plot.data$y, na.rm=T), reversed=TRUE)
  p$series(name = 'Heatmap',
           data =  toJSONArray2(plot.data, json=FALSE),
           color = "#cccccc",
           tooltip = list(headerFormat= '', pointFormat = paste('{point.Sample}<br/>{point.Accession}<br/>', 
                                                                dataType, ': {point.value}<br/>',
                                                                'Column: {point.x}<br/>Row: {point.y}', sep="")),
           dataLabels = list(
             enabled = showLabels,
             color = 'black',
             style = list(
               textShadow = 'none',
               HcTextStroke = NULL
             )
           ))
  
  # colorAxis is required for heat maps        
  p$addParams(colorAxis = 
                list(min = 0,
                     stops = list(
                                list(0.2, '#c4463a'),
                                list(0.6, '#fffbbc'),
                                list(1.2, '#3060cf')
                     )
                )
  )
  
  p$legend(align='right',
           layout='vertical',
           margin=0,
           verticalAlign='top',
           y=25,
           symbolHeight=320)
  
  if(!useWithShiny){
  #need to link additional JavaScript files here that do not ship with rCharts
  p$addAssets(js = 
                c("https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js",
                  "https://code.highcharts.com/highcharts.js",
                  "https://code.highcharts.com/highcharts-more.js",
                  "https://code.highcharts.com/modules/exporting.js",
                  "https://code.highcharts.com/modules/heatmap.js"
                )
  )
  }

  
  # save to standalone HTML page
  #p$save(destfile = 'heatmap.html')
  
  return(p)
}