library(rCharts)
library(rjson)
library(dplyr)

highcharts.scatterplot.plate <- function(plot.data, show.error=F){
        
    colnames(plot.data)[1:4] <- c("x", "y", "sem", "group")
        
    # create Highcharts object
    p <- Highcharts$new()
    
    # use type='heatmap' for heat maps
    p$chart(zoomType = "x")    
    p$xAxis(min = min(plot.data$x), max = max(plot.data$x))
    p$yAxis(min = min(plot.data$y), max = max(plot.data$y))
    
    #split data by group
    plot.data.split <- split(plot.data, plot.data$group)
    series <- foreach(group = names(plot.data.split), .combine=append) %do% {
        group.data <- plot.data.split[[group]]
        if(nrow(group.data) == 0) return(NULL)
        list(
            list(name = group,
                 type = 'scatter',
                 data =  toJSONArray2(group.data, json=FALSE),             
                 tooltip = list(pointFormat = '{point.Sample} {point.Accession}')
            ),
            list(name = group,#paste(group, 'error', sep=" "),
                 type = 'errorbar',
                 data = foreach(i = 1:nrow(group.data)) %do% 
                    {                        
                        return(c(group.data$x[i],c(group.data$y[i] - group.data$sem[i],
                                 group.data$y[i] + group.data$sem[i])))
                    },
                 enableMouseTracking = FALSE
            )
        )
    }
    
    p$series(
        series
    )        
            
    # save to standalone HTML page
    #p$save(destfile = 'scatter.html')
    
    return(p)
}