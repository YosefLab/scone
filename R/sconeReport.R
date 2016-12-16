#' SCONE Report Browser: Browse Evaluation of Normalization Performance
#' 
#' This function opens a shiny application session for visualizing performance 
#' of a variety of normalization schemes.
#' 
#' @param x a \code{SconeExperiment} object
#' @param methods character specifying the normalizations to report.
#' @param qc matrix. QC metrics to be used for QC evaluation report. Required.
#' @param bio factor. A biological condition (variation to be preserved).
#'   Default NULL.
#' @param batch factor. A known batch variable (variation to be removed).
#'   Default NULL.
#' @param negcon character. Genes to be used as negative controls for 
#'   evaluation. These genes should be expected not to change according to the 
#'   biological phenomenon of interest. Default empty character.
#' @param poscon character. Genes to be used as positive controls for 
#'   evaluation. These genes should be expected to change according to the 
#'   biological phenomenon of interest. Default empty character.
#'   
#' @importFrom shiny fluidPage downloadLink br titlePanel tableOutput helpText
#'   downloadHandler column h5 h6 updateSelectInput plotOutput tabsetPanel
#'   selectInput renderTable p sidebarLayout sidebarPanel shinyApp mainPanel
#'   observeEvent renderPlot reactive fluidRow sliderInput tabPanel
#' @import visNetwork
#' @importFrom DT dataTableProxy dataTableOutput renderDataTable selectRows
#'   formatSignif
#' @import plotly
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom NMF aheatmap
#' @importFrom ggplot2 labs theme geom_point ylim ggplot geom_violin
#'   element_blank aes geom_bar coord_cartesian scale_fill_manual guides
#' @export
#' 
#' @return An object that represents the SCONE report app.
#'   
#' @examples
#' set.seed(101)
#' mat <- matrix(rpois(1000, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat)
#' res <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
#'            evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2, 
#'            bpparam = BiocParallel::SerialParam())
#' qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
#' rownames(qc) = colnames(mat)
#' colnames(qc) = c("NCOUNTS","NGENES")
#' \dontrun{
#' sconeReport(res,rownames(get_params(res)), qc = qc)
#' }
#' 
sconeReport = function(x, methods,
                       qc, 
                       bio = NULL, batch = NULL,
                       poscon = character(), negcon = character()){
  
  
  scone_res = list()
  
  # List of normalized matrices
  scone_res$normalized_data = lapply(as.list(methods), 
                                     FUN = function(z){
                                       get_normalized(x,method = z,log=TRUE)
                                     })
  names(scone_res$normalized_data) = methods
  
  # Parameter matrix
  scone_res$params = get_params(x)[methods,]
  
  # Merged score matrix
  scone_res$scores = cbind(get_scores(x),get_score_ranks(x))[methods,]
  colnames(scone_res$scores)[ncol(scone_res$scores)] = "mean_score_rank"
  
  
  ## ----- If NULL classifications, Replace with NA ------
  
  if(is.null(bio)){
    bio = factor(rep("NA",ncol(scone_res$normalized_data[[1]])))
  }
  
  if(is.null(batch)){
    batch = factor(rep("NA",ncol(scone_res$normalized_data[[1]])))
  }
  
  ## ----- Network Data for Visualization -----
  
  # Matrix nodes in scone_res
  leaves = rownames(scone_res$params) 
  
  # Parents
  lay1 = unique(gsub(",[^,]*$","",leaves))
  lay2 = unique(gsub(",[^,]*$","",lay1))
  lay3 = unique(gsub(",[^,]*$","",lay2))
  lay4 = unique(gsub(",[^,]*$","",lay3))
  all_nodes = c(leaves,lay1,lay2,lay3,lay4)
  rm(list = c("lay1","lay2","lay3","lay4"))
  
  # Collapse no-op names for all node names
  all_nodes = unique(gsub("^,*|,*$","",
                          gsub("\\,+",",",
                               gsub("none|no_uv|no_bio|no_batch","",
                                    all_nodes))))
  
  # Order by decreasing number of operations (root at the end)
  all_nodes = all_nodes[order(paste0(gsub(".+",",",all_nodes),
                                     gsub("[^,]*","",all_nodes)),
                              decreasing = TRUE)]
  
  # Collapse no-op names for all leaves
  leaves = gsub("^,*|,*$","",
                gsub("\\,+",",",
                     gsub("none|no_uv|no_bio|no_batch","",
                          leaves)))
  
  # Re-title nodes in scone_res and root (instead of collapsed names)
  titles = all_nodes
  names(titles) =  paste0("_",all_nodes)
  titles[paste0("_",leaves)] = rownames(scone_res$params)
  titles[titles == ""] = "none"
  
  # Create Nodes
  nodes <- data.frame(id = 1:length(all_nodes),
                      title = titles,
                      group = c("NoData",
                                "Loaded")[1 + 
                                            (all_nodes  %in% leaves)])
  
  # Create Edges
  edges = data.frame()
  for(i in 1:length(all_nodes)){
    parents = as.vector(sapply(all_nodes,
                               function(x){
                                 grepl(paste0("^",x),all_nodes[i])
                               }))
    parents = parents & (all_nodes != all_nodes[i])
    if(any(parents)){
      edges = rbind(edges,cbind(min(which(parents )),i))
    }
  }
  colnames(edges)= c("from","to")
  
  # Re-map indices according to performance
  new_title_list = rownames(scone_res$params) # Ordered by perfomance
  new_title_list = c(new_title_list,
                     as.character(
                       nodes$title[!nodes$title %in%
                                     new_title_list])
  )
  old_ids = nodes[paste0("_",
                         gsub("^,*|,*$","",
                              gsub("\\,+",",",
                                   gsub("none|no_uv|no_uv|no_bio|no_batch",
                                        "",new_title_list)))),]$id
  old_nodes = nodes
  old_edges = edges
  for(i in 1:length(old_ids)){
    nodes$id[old_nodes$id == old_ids[i]] = i
    edges[old_edges == old_ids[i]] = i
  }
  nodes = nodes[order(nodes$id),]
  
  ## ----- Color Variable -----
  cc <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),brewer.pal(12,"Set3"),
          brewer.pal(8, "Accent"))
  
  ## ----- UI Definition -----
  ui <- fluidPage(
    titlePanel("SCONE Report Browser"),
    
    sidebarLayout(
      sidebarPanel(
        selectInput("norm_code", label = "Normalization",
                    choices = structure(as.list(rownames(scone_res$params)),
                                        names = rownames(scone_res$params)),
                    selected = rownames(scone_res$params)[1]),
        selectInput("color_code", label = "Stratify Plots by",
                    choices = list("Batch"="batch",
                                   "Biological Condition"="bio", 
                                   "PAM Cluster"="clust"),
                    selected = "batch"),
        sliderInput("dim", label = "Reduced Dimension (PCA)",
                    min=2, max=min(length(batch)-1,10),
                    value = 3),
        sliderInput("k", label = "Number of Clusters (PAM)",
                    min=2, max=10,
                    value = 2),
        helpText(paste0("Please note that with many cells (>500),",
                        " some plots may take several seconds ",
                        "(sometimes minutes) to appear.")),
        downloadLink('downloadData', 'Download')
      ),
      
      mainPanel(tabsetPanel(type = "tabs",
                            tabPanel("Overview",
                                     br(),
                                     visNetworkOutput("norm_net",
                                                      width = "650px"),
                                     br(),
                                     DT::dataTableOutput('tbl_score')
                            ),
                            tabPanel("PCA",
                                     br(),
                                     p(paste0("This panel shows principal ",
                                              "component analysis ",
                                              "(PCA) results",
                                              " on different subsets ",
                                              "of genes.")),
                                     br(),
                                     h6("Variance Explained (All Genes)"),
                                     p(paste0("Use this plot",
                                              " to decide on the ",
                                              "dimensionality of the reduced",
                                              " spaced used for evaluation.")),
                                     plotOutput("plot_scree",width = "650px",
                                                height = "400px"),
                                     br(),
                                     h6("2D (All Genes)"),
                                     plotOutput("plot_base",
                                                width = "650px",
                                                height = "450px"),
                                     br(),
                                     h6("3D Interactive (All Genes)"),
                                     plotlyOutput("plot3d_base",
                                                  width = "650px",
                                                  height = "650px"),
                                     br(),
                                     selectInput("gene_set", 
                                                 label = "Gene selection",
                                                 choices = 
                                                   list("Negative Controls"=
                                                          "neg",
                                                        "Positive Controls"=
                                                          "pos"),
                                                 selected = "neg"),
                                     h6("2D (Select Genes)"),
                                     p(paste0("Visualize cells by PCs",
                                              " of control gene sets.")),
                                     plotOutput("plot_select",
                                                width = "650px",
                                                height = "450px")
                            ),
                            tabPanel("QC",
                                     br(),
                                     p(paste0("This panel shows the absolute",
                                              " correlations between",
                                              " Principal",
                                              " Components (PCs) of",
                                              " expression",
                                              " matrix vs. PCs of the QC",
                                              " (Quality Control) measure ",
                                              "matrix. It also shows the PCA ",
                                              "of the QC measures. Note that",
                                              " the PCA score signs are not",
                                              " meaningful (e.g.,",
                                              " 'good' cells",
                                              " could have a low value of ",
                                              "PC1 and 'bad' cells a",
                                              " high value).")),
                                     br(),
                                     plotlyOutput('qccorPlot',
                                                  width = "650px",
                                                  height = "450px"),
                                     h5("PCA of QC metrics"),
                                     h6("2D"),
                                     plotOutput("plot_qc",
                                                width = "650px",
                                                height = "450px"),
                                     h6("3D Interactive"),
                                     plotlyOutput("plot3d_qc",
                                                  width = "650px",
                                                  height = "650px")
                            ),
                            tabPanel("Silhouette",
                                     br(),
                                     p(paste0("Silhouette Width per",
                                              " Sample and Contingency",
                                              " Tables")),
                                     br(),
                                     plotOutput("plotsil",
                                                width = "650px",
                                                height = "450px"),
                                     selectInput("cat1", 
                                                 label =
                                                   "Row Class",
                                                 choices =
                                                   list("Batch"=
                                                          "batch",
                                                        "Biological Condition"=
                                                          "bio", "PAM Cluster"=
                                                          "clust"),
                                                 selected = "batch"),
                                     selectInput("cat2",
                                                 label = "Column Class",
                                                 choices =
                                                   list("Batch"=
                                                          "batch",
                                                        "Biological Condition"=
                                                          "bio",
                                                        "PAM Cluster"=
                                                          "clust"),
                                                 selected = "bio"),
                                     tableOutput("cat_tab")
                            ),
                            tabPanel("Control Genes",
                                     br(),
                                     p(paste0("Heatmap of control genes, ",
                                              "colored by all",
                                              " three categories.")),
                                     br(),
                                     p("Positive controls:"),
                                     plotOutput("hmposcon",
                                                width = "650px",
                                                height = "450px"),
                                     p("Negative controls:"),
                                     plotOutput("hmnegcon",
                                                width = "650px",
                                                height = "450px")
                            ),
                            
                            tabPanel("Stratified PCA",
                                     br(),
                                     p(paste0("Violin plots showing",
                                              " conditional ",
                                              "distributions of",
                                              " selected Principal",
                                              " Component.")),
                                     br(),
                                     sliderInput("pcsel", 
                                                 label =
                                                   paste0("Selected",
                                                          " Principal ",
                                                          "Component"),
                                                 min=1, max=
                                                   min(length(batch)-1,10),
                                                 value = 1),
                                     h6("All Genes:"),
                                     plotOutput("violin_base",
                                                width = "650px",
                                                height = "450px"),
                                     selectInput("gene_set2", 
                                                 label = "Gene selection",
                                                 choices =
                                                   list("Negative Controls"=
                                                          "neg",
                                                        "Positive Controls"=
                                                          "pos"),
                                                 selected = "neg"),
                                     h6("Select Genes:"),
                                     plotOutput("violin_select",
                                                width = "650px",
                                                height = "450px")
                            ),
                            
                            tabPanel("Relative Log-Expression",
                                     br(),
                                     p(paste0("Relative ",
                                              "Log-Expression Plot ",
                                              "for top 100 Most ",
                                              "Variable Genes.")),
                                     plotOutput("rle",
                                                width = "850px",
                                                height = "450px")
                            )
      ))))
  
  server <- function(input, output, session) {
    
    ## ------ Overview Tab ------
    
    # Render Network
    output$norm_net <- renderVisNetwork({
      
      # Awk: Check if any non-loaded methods are included
      if(any(!all_nodes  %in% leaves)){
        visNetwork(nodes, edges,width = "100%",main = "Tree of Methods") %>% 
          visHierarchicalLayout(sortMethod = "directed" )  %>% 
          visGroups(groupname = "NoData", shape = "dot", 
                    size = 10, color = list(background = "lightgrey",
                                            border = "darkblue",
                                            highlight = 
                                              list(background =
                                                     "black",
                                                   border = "red")), 
                    shadow = list(enabled = TRUE)) %>% 
          visGroups(groupname = "Loaded", shape = "dot", 
                    size = 20, color = list(background = "lightblue",
                                            border = "darkblue", 
                                            highlight = 
                                              list(background = "cyan",
                                                   border = "red")),
                    shadow = list(enabled = TRUE)) %>%
          visEdges(shadow = TRUE, arrows = list(to = list(enabled = TRUE)),
                   color = list(color = "darkblue", highlight = "red")) %>%
          visOptions(nodesIdSelection =
                       list(enabled = TRUE, 
                            values = nodes$id[nodes$group == "Loaded"],
                            selected =
                              nodes$id[nodes$title ==
                                         rownames(scone_res$params)[1]])) %>%
          visLegend(width = 0.1, position = "right", main = "Status")
      }else{
        visNetwork(nodes, edges,width = "100%",main = "Tree of Methods") %>% 
          visHierarchicalLayout(sortMethod = "directed" )  %>% 
          visEdges(shadow = TRUE, arrows = list(to = list(enabled = TRUE)),
                   color = list(color = "darkblue", highlight = "red")) %>%
          visOptions(nodesIdSelection =
                       list(enabled = TRUE, 
                            values = nodes$id[nodes$group == "Loaded"],
                            selected = 
                              nodes$id[nodes$title ==
                                         rownames(scone_res$params)[1]]))
      }
      
    })
    
    # Render Table
    output$tbl_score <- DT::renderDataTable({
      DT::datatable(scone_res$scores,
                    extensions = 'Buttons',
                    selection = list( mode = 'single',
                                      target = 'row',
                                      selected = c(1)),
                    options=
                      list(columnDefs =
                             list(list(visible=FALSE,
                                       targets=1:(ncol(scone_res$scores)-1))),
                           dom = 'Bfrtip', 
                           buttons = list(
                             list(extend = 'colvis',
                                  columns = c(1:(ncol(scone_res$scores)-1))))
                      )) %>%
        DT::formatSignif(columns=c(1:(ncol(scone_res$scores))), digits=3)
    })
    
    # Update Menu Upon Network Selection of (Loaded) Normalization
    observeEvent(input$norm_net_selected,{
      if(is.character(input$norm_net_selected)){ # probably unnecessary 
        if(!input$norm_net_selected == ""){ # probably unnecessary
          if(as.integer(input$norm_net_selected) %in%
             nodes$id[nodes$group == "Loaded"]){
            updateSelectInput(
              session,
              "norm_code",
              selected = as.character(
                nodes$title[nodes$id ==
                              as.integer(input$norm_net_selected)]
              ))
          }
        }
      }
    })
    
    # Update Menu Upon Table Selection
    observeEvent(input$tbl_score_rows_selected,{
      if(length(input$tbl_score_rows_selected) > 0 ){ # probably unnecessary 
        updateSelectInput(session,
                          "norm_code", 
                          selected = as.character(
                            nodes$title[nodes$id ==
                                          as.integer(
                                            input$tbl_score_rows_selected
                                          )]
                          )
        )
      }
    })
    
    # Upon Menu Selection of Normalization
    observeEvent(input$norm_code,{
      
      # Update Network 
      if(is.character(input$norm_net_selected)){
        
        # If selection is empty, then no valid node is selected.
        if(input$norm_net_selected == ""){
          visNetworkProxy("norm_net") %>%
            visSelectNodes(id =
                             array(nodes$id[nodes$title == input$norm_code]))
        }else{
          
          # If selection is different from menu, then network must be updated.
          if(as.character(nodes$title[nodes$id == 
                                      as.integer(
                                        input$norm_net_selected
                                        )]) != input$norm_code){
            visNetworkProxy("norm_net") %>%
              visSelectNodes(id =
                               array(nodes$id[nodes$title == input$norm_code]))
          }
        }
      }
      
      # Update Table
      
      # If selection is empty, then no row is selected.
      if(length(input$tbl_score_rows_selected) == 0){
        dataTableProxy("tbl_score") %>%
          selectRows(nodes$id[nodes$title == input$norm_code])
      }else{
        
        # If selection is different from menu, then row must be updated.
        if(input$tbl_score_rows_selected !=
           nodes$id[nodes$title == input$norm_code]){
          dataTableProxy("tbl_score") %>%
            selectRows(nodes$id[nodes$title == input$norm_code])
        }
      }
      
    })
    
    ## ------ PCA Tab ------
    
    normle <- reactive({
      as.matrix(scone_res$normalized_data[[input$norm_code]])
    })
    
    strat_col <- reactive({
      switch(input$color_code,
             bio = bio,
             batch = batch,
             clust = pam_obj()$clust)
    })
    
    pc_col <- reactive({
      cc[strat_col()]
    })
    
    pc_gene <- reactive({
      switch(input$gene_set,
             pos = intersect(poscon,rownames(normle())),
             neg = intersect(negcon,rownames(normle())))
    })
    
    pc_obj_base <- reactive({
      prcomp(t(normle()),center = TRUE, scale = TRUE)
    })
    
    pc_obj_select <- reactive({
      if(length(pc_gene()) > 0){
        prcomp(t(normle()[pc_gene(),]),center = TRUE, scale = TRUE)
      }else{
        list()
      }
    })
    
    pam_obj <- reactive({
      cluster::pam(x = pc_obj_base()$x[,1:input$dim],k = input$k)
    })
    
    output$plot_scree <- renderPlot({
      plot(pc_obj_base()$sdev[1:min(input$dim*2,
                                    length(batch)-1)]^2/
             sum(pc_obj_base()$sdev^2), 
           typ = "l",
           ylab = "Proportion of Variance", 
           xlab = "PC Index",
           log = "y")
      abline(v = input$dim, col = "red",lty = 2)
    })
    
    output$plot_base <- renderPlot({
      par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
      plot(pc_obj_base()$x[,1:2],col =  pc_col(), pch = 16)
      legend("topright", inset=c(-0.3,0), 
             legend=levels(factor(strat_col())), 
             fill = cc[sort(unique(factor(strat_col())))])
    })
    
    output$plot3d_base <- renderPlotly({
      PC1 <- PC2 <- PC3 <- NULL
      df <- setNames(data.frame(pc_obj_base()$x[,1:3]),
                     c("PC1", "PC2", "PC3"))
      plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3,
              type = "scatter3d", mode = "markers",
              marker = list(color=pc_col() ))
    })
    
    output$plot_select <- renderPlot({
      
      par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
      
      if(length(pc_gene()) > 0){
        
        plot(pc_obj_select()$x[,1:2],
             col =  pc_col(), pch = 16)
        legend("topright", inset=c(-0.3,0), 
               legend=levels(factor(strat_col())),
               fill = cc[sort(unique(factor(strat_col())))])
        
      }else{
        
        plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
        text(0,labels = "Control gene set is empty.")
        
      }
      
    })
    
    pam_obj <- reactive({
      cluster::pam(x = pc_obj_base()$x[,1:input$dim],k = input$k)
    })
    
    ## ------ QC Tab ------
    
    cor_qc <- reactive({
      abs(cor(qc,pc_obj_base()$x))
    })
    
    output$qccorPlot <- renderPlotly({
      metric <- value <- PC <- NULL
      
      df = melt(cor_qc()[,1:input$dim])
      colnames(df) = c("metric","PC","value")
      df$metric = factor(df$metric,levels = colnames(qc))
      df$PC = factor(df$PC,
                     levels = paste0("PC",1:input$dim))
      p <- ggplot(data = df, aes(x = PC, 
                                 fill = metric, 
                                 weight=value)) +
        geom_bar(position = "dodge") + 
        ylim(0, 1) + 
        labs(x="PC of Expression",
             y="Absolute Correlation") + 
        theme(legend.title=element_blank())
      ggplotly(p)
    })
    
    pc_obj_qc <- reactive({
      prcomp(as.matrix(qc),center = TRUE, scale = TRUE)
    })
    
    output$plot_qc <- renderPlot({
      par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
      plot(pc_obj_qc()$x[,1:2],col =  pc_col(), pch = 16)
      legend("topright", inset=c(-0.3,0),
             legend=levels(factor(strat_col())),
             fill = cc[sort(unique(factor(strat_col())))])
    })
    
    output$plot3d_qc <- renderPlotly({
      if(ncol(pc_obj_qc()$x) >= 3){
        PC1 <- PC2 <- PC3 <- NULL
        df <- setNames(data.frame(pc_obj_qc()$x[,1:3]),
                       c("PC1", "PC2", "PC3"))
        plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, 
                type = "scatter3d", 
                mode = "markers",
                marker = list(color=pc_col() ))
      }else{
        plot_ly(data.frame(), type = "scatter3d", mode = "markers")
      }
    })
    
    ## ------ Silhouette Tab ------
    
    sil_obj <- reactive({
      cluster::silhouette(x = as.numeric(strat_col()),
                          dist = dist(pc_obj_base()$x[,1:input$dim]))
    })
    
    output$plotsil <- renderPlot({
      
      par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
      
      if(length(unique(strat_col())) > 1){
        
        sil = sil_obj()
        o1 = order(sil[,3])
        o = o1[order(-sil[,1][o1])]
        barplot(sil[,3][o], main = paste0("ASW = ",
                                          signif(mean(sil[,3]),3)),
                xlab = "silhouette width",horiz=TRUE,
                xlim = 1.5*range(c(sil[,3],-sil[,3])),
                col = pc_col()[o],
                border = pc_col()[o])
        legend("topright", inset=c(-0.3,0), 
               legend=levels(factor(strat_col())), 
               fill = cc[sort(unique(factor(strat_col())))])
        
      }else{
        
        plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
        text(0,labels = "Stratify plots by a multi-level classification.")
        
      }
      
    })
    
    cat1 <- reactive({
      switch(input$cat1,
             bio = bio,
             batch = batch,
             clust = pam_obj()$clust)
    })
    
    cat2 <- reactive({
      switch(input$cat2,
             bio = bio,
             batch = batch,
             clust = pam_obj()$clust)
    })
    
    output$cat_tab <- renderTable({
      table(cat1(), cat2())      
    })
    
    
    ## ------ Control Genes Tab ------
    
    silo <- reactive({
      sil = sil_obj()
      o1 = order(sil[,3])
      o1[order(-sil[,1][o1])]    
    })
    
    output$hmnegcon <- renderPlot({
      
      if(length(negcon) > 0){
        
        if(length(unique(strat_col())) > 1){
          
          aheatmap(normle()[negcon,], Colv = silo(), 
                   annCol = list(batch = batch, bio = bio,
                                 pam =  as.factor(pam_obj()$clust)),
                   annColors = list(batch = cc, bio = cc, pam =  cc))
          
        }else{
          
          plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
          text(0,labels = "Stratify plots by a multi-level classification.")
          
        }
        
      }else{
        
        plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
        text(0,labels = "Control gene set is empty.")
        
      }
      
    })
    
    output$hmposcon <- renderPlot({
      
      if(length(poscon) > 0){
        
        if(length(unique(strat_col())) > 1){
          
          aheatmap(normle()[poscon,], Colv = silo(), 
                   annCol = list(batch = batch, 
                                 bio = bio, 
                                 pam = as.factor(pam_obj()$clust)),
                   annColors = list(batch = cc, bio = cc, pam =  cc))
          
        }else{
          
          plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
          text(0,labels = "Stratify plots by a multi-level classification.")
          
        }
        
      }else{
        
        plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
        text(0,labels = "Control gene set is empty.")
        
      }
      
    })
    
    ## ------ Stratified PCA Tab ------
    
    output$violin_base <- renderPlot({
      
      Class = factor(strat_col())
      Val = pc_obj_base()$x[,input$pcsel] 
      
      ggplot(data.frame(Class,Val ),aes(x = Class,y = Val))   +  
        geom_violin(scale = "width", trim = TRUE, aes(fill = Class))+ 
        labs(x = "Plot Stratification", y = "PC") +
        coord_cartesian(ylim = max(abs(range(Val)))*c(-1.5,1.5))  +
        scale_fill_manual(values=cc[sort(unique(factor(strat_col())))]) +
        geom_point(colour = "black") + guides(fill=FALSE)
      
    })
    
    pc_gene2 <- reactive({
      switch(input$gene_set2,
             pos = intersect(poscon,rownames(normle())),
             neg = intersect(negcon,rownames(normle())))
    })
    
    pc_obj_select2 <- reactive({
      
      if(length(pc_gene2()) > 0){
        
        prcomp(t(normle()[pc_gene2(),]),center = TRUE, scale = TRUE)
        
      }else{
        list()
      }
      
    })
    
    output$violin_select <- renderPlot({
      
      if(length(pc_gene2()) > 0){
        
        Class = factor(strat_col())
        Val = pc_obj_select2()$x[,input$pcsel] 
        
        ggplot(data.frame(Class,Val ),aes(x = Class,y = Val))   +  
          geom_violin(scale = "width", trim = TRUE, aes(fill = Class))+ 
          labs(x = "Plot Stratification", y = "PC") +
          coord_cartesian(ylim = max(abs(range(Val)))*c(-1.5,1.5))  +
          scale_fill_manual(values=cc[sort(unique(factor(strat_col())))]) +
          geom_point(colour = "black") + guides(fill=FALSE)
        
      }else{
        
        plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
        text(0,labels = "Control gene set is empty.")
        
      }
    })
    
    ## ------ Relative Log-Expression ------
    
    output$rle <- renderPlot({
      
      par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
      
      if(length(unique(strat_col())) > 1){
        
        vars = apply(exp(normle()),1,var)
        is_var = rank(-vars) <= 100
        median <- apply(normle()[is_var,], 1, median)
        rle <- apply(normle()[is_var,], 2, function(x) x - median)
        boxplot(rle[,rev(silo())],col = pc_col()[rev(silo())], 
                outline = FALSE, names = rep("",ncol(rle)))
        abline(h=0, lty=2)
        legend("topright", inset=c(-0.2,0),
               legend=levels(factor(strat_col())), 
               fill = cc[sort(unique(factor(strat_col())))])
        
      }else{
        
        plot(0, type = 'n',xlab = "",ylab = "",xaxt = 'n',yaxt = 'n')
        text(0,labels = "Stratify plots by a multi-level classification.")
        
      }
      
    })
    
    ## ----- Download Button -----
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('scone_out-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        datt = cbind(as.character(bio),
                     as.character(batch),
                     as.character(pam_obj()$clust))
        colnames(datt) = c("Class-Bio","Class-Batch","Class-Pam")
        rownames(datt) = colnames(normle())
        datt = rbind(t(datt),normle())
        write.csv(datt, con)
      }
    )
    
  }
  
  # Shiny App
  shinyApp(ui = ui, server = server)
  
}
