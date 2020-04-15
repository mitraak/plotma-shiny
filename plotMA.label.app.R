library(ggplot2)
library(ggrepel)
library(DT)
library(dplyr)

# Increase max upload size to 30 MB
options(shiny.maxRequestSize = 30*1024^2)

# Plot a MA plot labeled with selected genes
#
# @param res.list data.frame (or list) with log2FoldChange and padj columns (or elements).
#        Row names of the data.frame should be gene IDs/symbols and match the format of lab.genes
# @param fdr.thres FDR threshold for defining statistical significance
# @param fc.thres log2FoldChange cutoff for defining statistical significance
# @param fc.lim User-defined limits for y-axis (log2FC). If NULL, this is defined as
#        the (floor, ceiling) of range(res$log2FoldChange)
# @param genes.to.label Genes to label on the MA plot. NULL, by default
# @param col Column of `res` in which to look for `genes.to.label`. If NULL,
#        rownames are used.
#
# @return Handle to ggplot
plotMA.label <- function(res,
                         fdr.thres=0.1,
                         fc.thres=0,
                         fc.lim=NULL,
                         genes.to.label=NULL
                         ){
  genes.to.label <- as.character(genes.to.label)
  nna <- sum(is.na(genes.to.label))
  if (nna > 0){
      warning(paste("Removing", nna, "NAs from gene list"))
      genes.to.label <- genes.to.label[!is.na(genes.to.label)]
  }
  # convert res to data frame
  res <- data.frame(res)

  # if y limits not specified
  if(is.null(fc.lim)){
    fc.lim <- range(res$log2FoldChange, na.rm=TRUE)
    fc.lim[1] <- floor(fc.lim[1])
    fc.lim[2] <- ceiling(fc.lim[2])
  }

  # change NA symbols to gene ID
  if('symbol' %in% colnames(res)){
    df <- res %>%
      mutate(geneid=rownames(res)) %>%
      mutate(symbol=as.character(symbol)) %>%
      mutate(symbol = replace(symbol, is.na(symbol), geneid[is.na(symbol)]))
  } else {
    df$symbol <- rownames(res)
  }
  # create column with plotting character based on fc.lim
  # change plotted values for those outside plot limits
  # to the limits
  df$shape <- 'in'
  df <- df %>%
    mutate(log2FoldChange = replace(log2FoldChange, is.na(log2FoldChange), 0)) %>%
    mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange > fc.lim[2], fc.lim[2])) %>%
    mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange < fc.lim[1], fc.lim[1])) %>%
    mutate(shape = replace(shape, log2FoldChange == fc.lim[2], 'above')) %>%
    mutate(shape = replace(shape, log2FoldChange == fc.lim[1], 'below')) %>%
    mutate(shape = as.factor(shape))

  # change colors for DE genes
  df$significant <- 'no'
  df <- df %>%
    mutate(significant = replace(significant,
                                 padj < fdr.thres & !is.na(padj) & abs(log2FoldChange) >= fc.thres,
                                 'yes')) %>%
    mutate(significant=as.factor(significant))

  p <- ggplot(df, aes(baseMean, log2FoldChange, color=significant, shape=shape)) +
    geom_point(alpha=0.8) + scale_x_log10()
  p <- p + scale_color_manual(breaks=c('no','yes'),
                              values=c('gray40','red'), guide='legend')
  p <- p + scale_shape_manual(breaks=c('in','above','below'),
                              values=c(16, 2, 6), guide=FALSE)
  p <- p + theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  p <- p + geom_hline(yintercept = 0, col="red", size=2, alpha=0.5)


  if(!is.null(genes.to.label)){
    # data frame with labels
    label.list <- df %>% filter(symbol %in% genes.to.label)

    # add labels
    p <- p + geom_point(data=label.list, col="black", pch=1, size=3)
    p <- p + geom_label_repel(data=label.list, aes(label=symbol),
                              fontface="italic", size=6, show.legend=FALSE)
  }
  return(p)

}

ui <- fluidPage(
                # take parameters that will be passed to plotMA.label
                # using Shiny

                # side bar
                sidebarPanel(
                             numericInput("fdr.thres", label = h3("FDR threshold"), value = 0.01),
                             numericInput("fc.thres", label = h3("log2FoldChange threshold"), value = 0),
                             sliderInput("fc.lim", label=h3("log2FoldChange limits"),
                                         min=-20, max=20, value=c(-5,5)),

                             textInput("genes.to.label",
                                       label = h3("Genes to label (comma-separated symbols)"),
                                       value = "SPP1, SMAD3, NFKB1"),

                             h3("Genes from click"),
                             tags$style(type='text/css', '#clickgenes {background-color:white; color:red; border:2px solid white;}'),
                             textOutput("clickgenes"), br(),
                             actionButton("reset.clickgenes", label='Reset')),

               # main panel
               mainPanel(
                         # Choose file with DE results
                         fileInput('res.file', label=h2('Upload file')),

                         plotOutput('maplot', height='800px', brush='plot_brush', click='plot_click'),
                         downloadButton('maplot_download', label = 'Download plot'),
                         fluidRow(
                                  column(10, dataTableOutput('table'))
                                  ))
)

server <- function(input, output, session){
    genes.to.label <- reactive({
      toupper(unlist(strsplit(input$genes.to.label, '\\,\\s*')))
    })
    
    # reactive element that returns table
    res.raw <- reactive({
      if(is.null(input$res.file)) return(NULL)
      return(read.table(input$res.file$datapath, sep='\t',
                        header=TRUE, row.names=1, quote=''))
    })
    
    # reactiveValue
    # initialized to point outside plot limits
    coords <- reactiveValues(x=0, y=0)
    
    # observe plot click and update accordingly
    observe({
      input$plot_click
      isolate({
        coords$x <- c(coords$x, input$plot_click$x)
        coords$y <- c(coords$y, input$plot_click$y)
      })
    })
  
    # get gene name(s) nearest to plot click
    genes.from.click <- eventReactive({
      input$plot_click
      input$reset.clickgenes
    }, {
      if(!is.null(input$plot_brush)) return(NULL)
      res <- res.raw()
      genes <- NULL
      for(i in 1:length(coords$x)){
        click.x <- coords$x[i]
        click.y <- coords$y[i]
        
        # adaptive xdel since it is on a log10 scale
        if(click.x != 1){
          xdel <- 0.1*(10^(ceiling(log10(click.x))))
        } else {
          xdel <- 1
        }
        
        # ydel is 2% of visible range of log2FoldChange
        ydel <- 0.02*min(diff(input$fc.lim),
                         diff(range(res$log2FoldChange, na.rm=TRUE)))
        
        # get filtered indices in area: click coordinates with wiggle room
        idx.x <- res$baseMean > (click.x - xdel) & res$baseMean < (click.x + xdel)
        idx.y <- res$log2FoldChange > (click.y - ydel) & res$log2FoldChange < (click.y + ydel)
        
        # labels are symbols. If symbol = NA, then label=rowname
        labels <- as.character(res$symbol)
        labels[is.na(labels)] <- rownames(res)[is.na(labels)]
        
        # get genes within area of interest
        idx <- idx.x & idx.y
        g <- labels[idx]
        
        # find closest point, return first for ties
        if(length(g) > 1){
          df.subset <- res[idx,]
          x.dist <- df.subset$baseMean - click.x
          y.dist <- df.subset$log2FoldChange - click.y
          total.dist <- x.dist^2 + y.dist^2
          min.idx <- which(total.dist == min(total.dist))
          
          genes <- c(genes, g[min.idx][1])
        } else {
          genes <- c(genes, g)
        }
      }
      return(genes)
    }, ignoreNULL = FALSE)
    
    # maintain list of genes from click
    output$clickgenes <- renderText({
      g <- genes.from.click()
      if(sum(!is.na(g)) == 0){
        return(' ')
      } else {
        return(paste(g[!is.na(g)], collapse=', '))
      }
    })
    
    # reset click genes
    observeEvent(input$reset.clickgenes, {
      coords$x <- 0
      coords$y <- 0
    })
    
    # reactive element to generate MA plot
    maplot <- reactive({
      if(!is.null(input$res.file)){
        res <- res.raw()
      }
      p <- plotMA.label(res,
                        fdr.thres=input$fdr.thres,
                        fc.thres=input$fc.thres,
                        fc.lim=input$fc.lim,
                        genes.to.label=c(genes.to.label(),
                                         genes.from.click()))
      p + theme(text=element_text(size=20))
    })
    
    # element to output MA plot
    output$maplot <- renderPlot({
      if(is.null(input$res.file)) return(NULL)
      maplot()
      
      # zoom into brushed area
      #if(!is.null(input$plot_brush)){
      #  coords <- input$plot_brush
      #  p + xlim(coords$xmin, coords$xmax) +
      #    ylim(coords$ymin, coords$ymax) +
      #    scale_x_log10()
      #} else {
      #  p
      #}
    })
    
    # download button for plot
    output$maplot_download <- downloadHandler(
      filename = function(){
        paste0('maplot.pdf')
      },
      content = function(file){
        ggsave(file, maplot(), width=12, height=10)
      }
    )
    
    # output data table of selected area
    output$table <- renderDataTable({
      if(is.null(input$plot_brush)) return(NULL)
      res <- res.raw()
      lims <- input$plot_brush
      idx.x <- res$baseMean > lims$xmin & res$baseMean < lims$xmax
      idx.y <- res$log2FoldChange > lims$ymin & res$log2FoldChange < lims$ymax
      cbind(genes=rownames(res),res)[idx.x & idx.y,]
    })
}

shinyApp(ui = ui, server = server)
