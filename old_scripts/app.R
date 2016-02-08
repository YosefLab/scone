library(shiny)
library(EDASeq)
library(clusterCells)
load("../shiny/tmp.rda")

library(RColorBrewer)
cc <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),brewer.pal(12,"Set3"),
        brewer.pal(8, "Accent"))

choose_color <- function(input) {
  if(input$color_code=="type") {
    cols <- cc[type]
  } else {
    cols <- cc[batch]
  }
  return(cols)
}

server <- function(input, output) {
    output$boxplot <- renderPlot({
      set.seed(19883)
      if(input$downsample) {
        idx <- sample(1:100)
      } else {
        idx <- 1:ncol(norm)
      }      
      boxplot(log_norm[,idx], boxcol=cc[get(input$color_code)[idx]], medlwd=5, col=cc[get(input$color_code)[idx]],
              outline=FALSE, names=FALSE, ylab="log(count+1)", main="Log normalized counts")
      lines(apply(log_norm[,idx], 2, median))
    }
    )
    
    output$rle <- renderPlot({
      set.seed(19883)
      if(input$downsample & ncol(norm)>100) {
        idx <- sample(1:100)
      } else {
        idx <- 1:ncol(norm)
      }
      boxplot(rle[,idx], boxcol=cc[get(input$color_code)[idx]], medlwd=5, col=cc[get(input$color_code)[idx]],
              outline=FALSE, names=FALSE, ylab="Relative Log Expression", main="Relative Log Expression (RLE)")
      abline(h=0, lty=2)
    })
    
    output$pca <- renderPlot(
      plot(pca$x[,1:2], col=choose_color(input), pch=19, main="All genes")
    )
    
    output$pca_pc <- renderPlot(
      plot(pca_pc$x[,1:2], col=choose_color(input), pch=19, main="Positive controls")
    )
    
    output$pca_nc <- renderPlot(
      plot(pca_nc$x[,1:2], col=choose_color(input), pch=19, main="Negative controls")
    )

    output$pca_mv <- renderPlot(
      if(input$variance=="variance") {
        plot(pca_var$x[,1:2], col=choose_color(input), pch=19, main="Top 100 most variable genes")        
      } else {
        plot(pca_cv$x[,1:2], col=choose_color(input), pch=19, main="Top 100 most variable genes")        
      } 
    )
    
    output$qc_plot <- renderPlot({
      barplot(abs(cors), beside = TRUE, col=cc[1:ncol(qc)], border=cc[1:ncol(qc)], ylim=c(0, 1),
              space=c(0, 4), main="Absolute correlation with QC scores")
      legend("topright", colnames(qc), fill=cc, border=cc, cex=.5)  
    }
    )

    output$qc_pca <- renderPlot({
      plot(qc_pca$x[,1:2], pch=19, col=choose_color(input), main="Quality PCA")
    }
    )
    
    output$qc_corr <- renderTable(abs(cors[,1:5]))
    
    output$pam_cc <- renderPlot({
      dualHeatmap(heatData = sub[[input$k]], main="Heatmap of co-clustering, based on 100 subsamplings",
                  clusterSamples=TRUE, clusterVars=TRUE, dual=FALSE, 
                  annLegend=FALSE, annCol = data.frame(Batch=batch, Type=type),
                  annColor=list(Batch=cc, Type=cc))
    })
    
    output$pam_type <- renderTable({
      table(pam[[input$k-1]]$clustering, type)      
    })
    
    output$pam_batch <- renderTable({
      table(pam[[input$k-1]]$clustering, batch)      
    })
    
    output$pam_pca <- renderPlot({
      plot(pca$x[,1:2], pch=19, col=cc[pam[[input$k-1]]$clustering], main="Expression PCA")
      legend("topright", as.character(1:input$k), fill=cc, cex=.5)
    })
    
    output$pam_qc <- renderPlot({
      plot(qc_pca$x[,1:2], pch=19, col=cc[pam[[input$k-1]]$clustering], main="Quality PCA")
      legend("topright", as.character(1:input$k), fill=cc, cex=.5)
    })
    
    output$heatmap_pc <- renderPlot({
      dualHeatmap(heatData = t(norm), main="Positive Control genes",
                  clusterSamples=TRUE, clusterVars=TRUE, dual=FALSE, whVars=poscon,
                  annLegend=FALSE, annCol = data.frame(Batch=batch, Type=type),
                  annColor=list(Batch=cc, Type=cc))
    })
    
    output$heatmap_nc <- renderPlot({
      dualHeatmap(heatData = t(norm), main="Negative Control genes",
                  clusterSamples=TRUE, clusterVars=TRUE, dual=FALSE, whVars=negcon,
                  annLegend=FALSE, annCol = data.frame(Batch=batch, Type=type),
                  annColor=list(Batch=cc, Type=cc))
    })

    output$heatmap_mv <- renderPlot({
      if(input$variance=="variance") {
        genes <- mv_var
      } else {
        genes <- mv_cv
      } 
      
      dualHeatmap(heatData = t(norm), main="Top 100 most variable genes",
                  clusterSamples=TRUE, clusterVars=TRUE, dual=FALSE, whVars=genes,
                  annLegend=FALSE, annCol = data.frame(Batch=batch, Type=type),
                  annColor=list(Batch=cc, Type=cc))
    })

    output$heatmap_pc1 <- renderPlot({
      genes <- names(sort(pca$rotation[,1], decreasing = TRUE)[1:100])
      dualHeatmap(heatData = t(norm), main="Top 100 genes by PC1 loadings",
                  clusterSamples=TRUE, clusterVars=TRUE, dual=FALSE, whVars=genes,
                  annLegend=FALSE, annCol = data.frame(Batch=batch, Type=type),
                  annColor=list(Batch=cc, Type=cc))
    })
    
    output$box_pc1_pc <- renderPlot({
      boxplot(tapply(pca_pc$x[,1], get(input$color_code), function(x) x), las=2, col=cc, ylab="PC1 score", main="PCA of Positive Controls")
    })

    output$box_pc1_nc <- renderPlot({
      boxplot(tapply(pca_nc$x[,1], get(input$color_code), function(x) x), las=2, col=cc, ylab="PC1 score", main="PCA of Negative Controls")
    })
    
    output$box_pc1_qc <- renderPlot({
      boxplot(tapply(qc_pca$x[,1], get(input$color_code), function(x) x), las=2, col=cc, ylab="PC1 score", main="PCA of quality scores")
    })
    
    output$curve_dropout <- renderPlot({
      plot(logistic(fitted[order(fnr$Alpha[1,]),1])~sort(fnr$Alpha[1,]), type='l', ylim=c(0, 1), col=rgb(0, 0, 0, 0.2), xlab="log average expression", ylab="Probability of drop-out", main="Drop-out Curves")
      for(i in 2:ncol(fitted)) {
        lines(logistic(fitted[order(fnr$Alpha[1,]),i])~sort(fnr$Alpha[1,]), type='l', col=rgb(0, 0, 0, 0.2))
      }
      lines(logistic(rowMeans(fitted)[order(fnr$Alpha[1,])])~sort(fnr$Alpha[1,]), col=2, lwd=2)
    })
        
    output$box_detection <- renderPlot({
      boxplot(tapply(colSums(norm>5), get(input$color_code), function(x) x), las=2, col=cc, ylab="Number of detected genes (> 5 reads)", main="Detected Genes")
    })

    output$box_dropout <- renderPlot({
      boxplot(tapply(colMeans(1-w), get(input$color_code), function(x) x), las=2, col=cc, ylab="Drop-out rate", main="Drop-out Rate")
    })
    
  }

ui <- fluidPage(
  titlePanel("In-depth Normalization Report"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("color_code", label = "Color code / stratify by",
                  choices = list("Batch"="batch", "Cell Type"="type"),
                  selected = "batch"),
      selectInput("variance", label = "Select most variable genes by",
                  choices = list("Variance"="variance", "Coefficient of variation"="cv"),
                  selected = "cv"),
      sliderInput("k", label = "Number of clusters",
                  min=2, max=10,
                  value = 7),
      checkboxInput("downsample", label="Downsample to 100 cells for boxplot panel", value=TRUE),
      helpText("RLE: relative log expression."),
      helpText("PCA: principal component analysis."),
      helpText("QC: quality control."),
      helpText("The 'Boxplot' panel shows the read count distribution and the RLE plot."),
      helpText("The 'PCA' panel shows the PCA of different subset of genes."),
      helpText("The 'QC' panel shows the absolute correlation between the first 10 principal components of the expression matrix and the QC measures. It also shows the PCA of the QC measures. Note that the PCA score signs are not meaningful (e.g., 'good' cells could have a low value of PC1 and 'bad' cells a high value)."),
      helpText("The 'Cluster stability' panel shows a heatmap of the PAM co-clustering based on 100 subsamples of the data. It shows the relation between the clusters obtained by PAM (on the original data) to cell type, batch, expression PCA, and QC PCA. Note that these results are reported for normalization evaluation only. A proper cluster analysis should be performed to find sub-populations in the data."),
      helpText("The 'Heatmaps' panel shows the log expression of the positive controls, negative controls, genes with highest PC1 loadings, and most variable genes."),
      helpText("The 'Confounding' panel shows the first PC of positive controls, negative controls, and QC stratified by either cell type or batch. This is useful to identify lower quality batches and potential confounding between quality and cell type."),
      helpText("The 'Drop-out' panel shows the number of detected genes and the drop-out rate stratified by either cell type or batch. It also shows the probability of drop-out as a function of the average expression"),
      helpText("The 'Number of clusters' slider controls the PAM clustering in the 'Cluster stability' tab."),
      helpText("Please note that with many cells (>500), some plots may take several seconds (sometimes minutes) to appear."),
      helpText("By default, we downsample the data to 100 cells for the boxplot visualization (only for the 'Boxplot' tab). Please uncheck it to see the full dataset (it may take some time).")
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  ## Add a "summary" tab (wait for Michael's text)
                  tabPanel("Boxplot",
                           plotOutput("boxplot"),
                           plotOutput("rle")),
                  tabPanel("PCA",
                           plotOutput("pca"),
                           plotOutput("pca_pc"),
                           plotOutput("pca_nc"),
                           plotOutput("pca_mv")),
                  tabPanel("QC",
                           plotOutput("qc_plot"),
                           tableOutput("qc_corr"),
                           plotOutput("qc_pca")),
                  tabPanel("Cluster stability",
                           plotOutput("pam_cc"),
                           fluidRow(
                             column(3,
                                    h3("Type"),
                                    br(),
                                    tableOutput("pam_type")),
                             column(3,
                                    h3("Batch"),
                                    br(),
                                    tableOutput("pam_batch"))),
                           plotOutput("pam_pca"),
                           plotOutput("pam_qc")),
                  tabPanel("Heatmaps",
                           plotOutput("heatmap_pc"),
                           plotOutput("heatmap_nc"),
                           plotOutput("heatmap_pc1"),
                           plotOutput("heatmap_mv")),
                  tabPanel("Confounding",
                           plotOutput("box_pc1_pc"),
                           plotOutput("box_pc1_nc"),
                           plotOutput("box_pc1_qc")),
                  tabPanel("Drop-out",
                           plotOutput("box_detection"),
                           plotOutput("box_dropout"),
                           plotOutput("curve_dropout"))
                  )
                  
    )))

shinyApp(ui = ui, server = server)
