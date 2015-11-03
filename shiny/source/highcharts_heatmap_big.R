library(rCharts)
library(rjson)

test <- subset(MCF7, Plate %in% c(1,2,3) & cellType=="BCSC")[,c("column", "row", "CTB", "Sample", "Replicate", "Plate")]
test$y <- (test$Replicate - 1) * 8 + test$row
test$x <- (test$Plate - 1) * 12 + test$column
test$value <- test$CTB
library(dplyr)
test <- arrange(test, x, y)

# create Highcharts object
p <- Highcharts$new()

# use type='heatmap' for heat maps
p$chart(zoomType = "xy", type = 'heatmap')
#p$credits(text = "Created with rCharts and Highcharts", href = "http://rcharts.io")
p$title(text='Sales per employee per weekday')
#p$xAxis(min = 1, max = 36, alternateGridColor="#aaaaaa", tickPositions=list(1,13,25), categories=list("Plate 1", "Plate 2"))
p$xAxis(min = 1, max = 36, plotBands = list( color = 'orange', from=1, to=12, label = list(text="Plate 1", y=-10)))
p$yAxis(min = 1, max = 24, reversed=TRUE, plotLines = list( color = 'green', value=8, label = list(text="Replicate 1")))
p$series(name = 'Heatmap',
         data =  toJSONArray2(test, json=FALSE),
         color = "#cccccc",
          tooltip = list(headerFormat= '', pointFormat = '{point.Sample}'))
#          dataLabels = list(
#            enabled = TRUE,
#            color = 'black',
#            style = list(
#              textShadow = 'none',
#              HcTextStroke = NULL
#            )
#          ))

# colorAxis is required for heat maps        
p$addParams(column=list(pointPlacement = -0.5), colorAxis = 
              list(min = 0,
                   stops = list(
                              list(0.2, '#c4463a'),
                              list(0.6, '#fffbbc'),
                              list(1.2, '#3060cf')
                   )
                   #minColor='#FFFFFF',
                   #maxColor='#7cb5ec'
              )
)
p$addParams(exporting = list(enabled = TRUE, width = 102))

p$legend(align='right',
         layout='vertical',
         margin=0,
         verticalAlign='top',
         y=25,
         symbolHeight=320)

# custom tooltip
#p$tooltip(formatter = "#! function() { return '<b>' + this.series.xAxis.categories[this.point.x] + '</b> sold <br><b>' +
                 #   this.point.value + '</b> items on <br><b>' + this.series.yAxis.categories[this.point.y] + '</b>'; } !#")

# need to link additional JavaScript files here that do not ship with rCharts
p$addAssets(js = 
              c("https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js",
                "https://code.highcharts.com/highcharts.js",
                "https://code.highcharts.com/highcharts-more.js",
                "https://code.highcharts.com/modules/exporting.js",
                "https://code.highcharts.com/modules/heatmap.js"
              )
)

# save to standalone HTML page
p$save(destfile = 'heatmap.html')

# open chart in browser
print(p)