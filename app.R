# Load libraries
library(shiny)
library(DT)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# ---- Helper for custom combined plot ----
plot_combined_enrichment <- function(up_df, down_df, title, out_file = NULL) {
  # --- Skip if both empty ---
  if ((is.null(up_df) || nrow(up_df) == 0) && (is.null(down_df) || nrow(down_df) == 0)) return(NULL)
  
  # --- Process Up ---
  if (!is.null(up_df) && nrow(up_df) > 0) {
    up_df <- head(up_df[order(up_df$p.adjust), ], 10)
    up_df$Term <- up_df$Description
    up_df$Direction <- "Up"
    up_df$logP <- -log10(up_df$p.adjust)
  } else {
    up_df <- data.frame(Term = character(), logP = numeric(), Direction = character())
  }
  
  # --- Process Down ---
  if (!is.null(down_df) && nrow(down_df) > 0) {
    down_df <- head(down_df[order(down_df$p.adjust), ], 10)
    down_df$Term <- down_df$Description
    down_df$Direction <- "Down"
    down_df$logP <- -(-log10(down_df$p.adjust)) # make negative
  } else {
    down_df <- data.frame(Term = character(), logP = numeric(), Direction = character())
  }
  
  df_combined <- rbind(up_df, down_df)
  if (nrow(df_combined) == 0) return(NULL)
  
  df_combined$Direction <- factor(df_combined$Direction, levels = c("Up", "Down"))
  df_combined <- df_combined[order(df_combined$Direction, -abs(df_combined$logP)), ]
  df_combined$Term <- factor(df_combined$Term, levels = rev(df_combined$Term))
  
  p <- ggplot(df_combined, aes(x = Term, y = logP, fill = Direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#E64B35FF", "Down" = "#4DBBD5FF")) +
    labs(
      title = title,
      y = expression(-log[10]("Adjusted P-value")),
      x = NULL,
      fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10, color = "black", face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold"),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    expand_limits(y = c(min(df_combined$logP) * 1.2, max(df_combined$logP) * 1.2))
  
  if (!is.null(out_file)) {
    ggsave(out_file, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
  }
  p
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Up/Down Gene Expression + GO Enrichment (Interactive & Downloadable)"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload file (3 columns: Gene | log2FC | Statistic):",
                accept = c(".txt", ".csv")),
      helpText("Columns must be: Gene | log2FC | p-value/FDR"),
      tags$hr(),
      
      sliderInput("fc_thr", "Log2 Fold Change threshold:", min = 0, max = 5, value = 0.58, step = 0.01),
      sliderInput("pval_thr", "Statistic threshold:", min = 0.0001, max = 0.2, value = 0.05, step = 0.001),
      tags$hr(),
      
      selectInput("organism", "Organism:",
                  choices = c("Human" = "human", "Mouse" = "mouse"), selected = "human"),
      selectInput("ont", "GO Ontology:",
                  choices = c("Biological Process (BP)" = "BP",
                              "Molecular Function (MF)" = "MF",
                              "Cellular Component (CC)" = "CC"),
                  selected = "BP"),
      selectInput("padj_method", "p-Value Adjustment Method:",
                  choices = c("BH", "BY", "holm", "bonferroni", "fdr"), selected = "BH"),
      sliderInput("p_cutoff", "p-Value Cutoff:", min = 0, max = 0.1, value = 0.05, step = 0.001),
      sliderInput("q_cutoff", "q-Value Cutoff:", min = 0, max = 0.5, value = 0.2, step = 0.01),
      numericInput("show_n", "Number of Categories to Show:", value = 10, min = 5, max = 30, step = 1),
      
      actionButton("filter_btn", "Apply thresholds & Run Analysis", class = "btn btn-primary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Up/Down Tables",
                 h4("Upregulated Genes"), DTOutput("up_table"),
                 h4("Downregulated Genes"), DTOutput("down_table")),
        
        tabPanel("GO Enrichment Plots",
                 h4("Upregulated GO Results"),
                 downloadButton("download_up_dot", "Download Dotplot (PNG)", class = "btn-info"),
                 plotOutput("up_go_dot", height = "500px"),
                 downloadButton("download_up_bar", "Download Barplot (PNG)", class = "btn-info"),
                 plotOutput("up_go_bar", height = "500px"),
                 tags$hr(),
                 h4("Downregulated GO Results"),
                 downloadButton("download_down_dot", "Download Dotplot (PNG)", class = "btn-info"),
                 plotOutput("down_go_dot", height = "500px"),
                 downloadButton("download_down_bar", "Download Barplot (PNG)", class = "btn-info"),
                 plotOutput("down_go_bar", height = "500px")),
        
        tabPanel("Combined Enrichment",
                 h4("Up + Down Combined Barplot"),
                 downloadButton("download_combined_bar", "Download Combined Barplot (PNG)", class = "btn-success"),
                 plotOutput("combined_bar", height = "600px"))
      )
    )
  )
)

# ---- SERVER ----
server <- function(input, output) {
  # --- Data reading ---
  dataInput <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    df <- switch(ext,
                 "txt" = read.delim(input$file$datapath, header = TRUE),
                 "csv" = read.csv(input$file$datapath, header = TRUE),
                 shiny::validate("Unsupported file type. Use .txt or .csv"))
    if (ncol(df) != 3)
      shiny::validate("File must have exactly 3 columns: Gene | log2FC | Statistic.")
    colnames(df) <- c("Gene", "log2FC", "Stat")
    df
  })
  
  filteredData <- eventReactive(input$filter_btn, {
    df <- dataInput()
    Up <- df[df$log2FC > input$fc_thr & df$Stat < input$pval_thr, ]
    Down <- df[df$log2FC < -input$fc_thr & df$Stat < input$pval_thr, ]
    list(Up = Up, Down = Down)
  })
  
  output$up_table <- renderDT({
    req(filteredData())
    datatable(filteredData()$Up, options = list(scrollX = TRUE))
  })
  output$down_table <- renderDT({
    req(filteredData())
    datatable(filteredData()$Down, options = list(scrollX = TRUE))
  })
  
  # --- Select organism ---
  getOrgDb <- reactive({
    if (input$organism == "human") org.Hs.eg.db else org.Mm.eg.db
  })
  
  # --- GO enrichment ---
  go_enrichment <- function(gene_list, orgdb, ont, p_cutoff, q_cutoff, padj_method) {
    if (nrow(gene_list) == 0) return(NULL)
    entrez <- bitr(gene_list$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
    if (is.null(entrez) || nrow(entrez) == 0) return(NULL)
    
    enrichGO(
      gene          = unique(entrez$ENTREZID),
      OrgDb         = orgdb,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = padj_method,
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff,
      readable      = TRUE
    )
  }
  
  go_results <- eventReactive(input$filter_btn, {
    orgdb <- getOrgDb()
    list(
      up = go_enrichment(filteredData()$Up, orgdb, input$ont, input$p_cutoff, input$q_cutoff, input$padj_method),
      down = go_enrichment(filteredData()$Down, orgdb, input$ont, input$p_cutoff, input$q_cutoff, input$padj_method)
    )
  })
  
  # --- DOTPLOTS ---
  output$up_go_dot <- renderPlot({
    res <- go_results()$up
    if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
    dotplot(res, showCategory = input$show_n, title = paste("Upregulated - GO", input$ont))
  })
  
  output$down_go_dot <- renderPlot({
    res <- go_results()$down
    if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
    dotplot(res, showCategory = input$show_n, title = paste("Downregulated - GO", input$ont))
  })
  
  # --- Custom BARPLOTS (-log10) ---
  custom_barplot <- function(enrich_res, title, n_show = 10) {
    df <- as.data.frame(enrich_res)
    if (nrow(df) == 0) return(NULL)
    
    df <- df %>% 
      arrange(p.adjust) %>%
      head(n_show) %>%
      mutate(negLog10 = -log10(p.adjust),
             Description = factor(Description, levels = rev(Description)))
    
    ggplot(df, aes(x = negLog10, y = Description)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal(base_size = 14) +
      labs(x = "-log10(adj. p-value)", y = NULL, title = title)
  }
  
  output$up_go_bar <- renderPlot({
    res <- go_results()$up
    if (is.null(res)) return(NULL)
    custom_barplot(res, paste("Upregulated - GO", input$ont), input$show_n)
  })
  
  output$down_go_bar <- renderPlot({
    res <- go_results()$down
    if (is.null(res)) return(NULL)
    custom_barplot(res, paste("Downregulated - GO", input$ont), input$show_n)
  })
  
  # --- Combined Up/Down plot ---
  output$combined_bar <- renderPlot({
    up_df <- as.data.frame(go_results()$up)
    down_df <- as.data.frame(go_results()$down)
    plot_combined_enrichment(up_df, down_df,
                             paste("Combined GO", input$ont))
  })
  
  # --- PNG Download Handlers ---
  plot_to_png <- function(plot_expr) {
    tmp <- tempfile(fileext = ".png")
    ggsave(tmp, plot = plot_expr, width = 8, height = 6, dpi = 300)
    tmp
  }
  
  output$download_up_dot <- downloadHandler(
    filename = function() paste0("Up_GO_", input$ont, "_dotplot.png"),
    content = function(file) {
      res <- go_results()$up
      if (is.null(res)) return(NULL)
      p <- dotplot(res, showCategory = input$show_n,
                   title = paste("Upregulated - GO", input$ont))
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  output$download_down_dot <- downloadHandler(
    filename = function() paste0("Down_GO_", input$ont, "_dotplot.png"),
    content = function(file) {
      res <- go_results()$down
      if (is.null(res)) return(NULL)
      p <- dotplot(res, showCategory = input$show_n,
                   title = paste("Downregulated - GO", input$ont))
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  output$download_up_bar <- downloadHandler(
    filename = function() paste0("Up_GO_", input$ont, "_barplot.png"),
    content = function(file) {
      res <- go_results()$up
      if (is.null(res)) return(NULL)
      p <- custom_barplot(res, paste("Upregulated - GO", input$ont))
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  output$download_down_bar <- downloadHandler(
    filename = function() paste0("Down_GO_", input$ont, "_barplot.png"),
    content = function(file) {
      res <- go_results()$down
      if (is.null(res)) return(NULL)
      p <- custom_barplot(res, paste("Downregulated - GO", input$ont))
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  output$download_combined_bar <- downloadHandler(
    filename = function() paste0("Combined_GO_", input$ont, "_barplot.png"),
    content = function(file) {
      up_df <- as.data.frame(go_results()$up)
      down_df <- as.data.frame(go_results()$down)
      plot_combined_enrichment(up_df, down_df,
                               paste("Combined GO", input$ont), out_file = file)
    }
  )
}

shinyApp(ui, server)
