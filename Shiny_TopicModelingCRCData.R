# Load required libraries

library(shiny)
library(shinydashboard)
library(curatedMetagenomicData)
library(tidyverse)
library(phyloseq)
library(LinDA)
library(topicmodels)
library(tidytext)
library(ldatuning)
library(cowplot)
library(viridis)
library(glue)
library(DT)            # For enhanced data table
library(shinycssloaders)  # For loading spinners

# Define UI for the dashboard
ui <- dashboardPage(
  dashboardHeader(title = "Microbiome Analysis Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "intro", icon = icon("info-circle")),
      menuItem("Data Overview", tabName = "data_overview", icon = icon("table")),
      menuItem("Topic Modeling", tabName = "topic_modeling", icon = icon("chart-bar")),
      menuItem("Differential Abundance", tabName = "diff_abundance", icon = icon("dna")),
      menuItem("Genera Analysis", tabName = "genera_analysis", icon = icon("bacteria")),
      menuItem("Selected Genera", tabName = "selected_genera", icon = icon("list")),
      menuItem("About", tabName = "about", icon = icon("question-circle"))
    )
  ),
  dashboardBody(
    tabItems(
      # Introduction tab
      tabItem(
        tabName = "intro",
        fluidRow(
          box(
            width = 12,
            h2("Microbiome Analysis Dashboard"),
            p("This dashboard provides interactive analysis of the microbiome data from the GuptaA_2019 study, including data visualization, topic modeling, and differential abundance analysis."),
            p("Please use the sidebar to navigate to different analysis sections."),
            br(),
            h3("Research Methods"),
            p("In this study, we compared the intestinal microbiota of patients with colorectal cancer (CRC) and healthy controls. The main methods include:"),
            tags$ul(
              tags$li(strong("Data preprocessing:"),"Filter and standardize the raw data."),
              tags$li(strong("Topic modeling:"), "Use the LDA algorithm to identify potential topic structures in the microbiome."),
              tags$li(strong("Differential abundance analysis:"), "Use the LinDA method to detect bacterial communities that are significantly different between CRC and the control group.")
            ),
            br(),
            p("Through these analyses, we hope to reveal the association between the microbiome and colorectal cancer and provide new insights into the diagnosis and treatment of the disease.")
          )
        )
      ),
      # Data Overview tab (Updated with DT::DTOutput)
      tabItem(
        tabName = "data_overview",
        fluidRow(
          box(
            width = 7,
            h2("Data Overview"),
            dataTableOutput("sample_table"),
            height = "600px"
          ),
          box(
            width = 5,
            title = "Case Status Distribution",
            withSpinner(plotOutput("case_status_plot", height = "480px")),
            height = "600px"
          )
        )
      ),
      # Topic Modeling tab
      tabItem(
        tabName = "topic_modeling",
        fluidRow(
          box(
            width = 6,
            height = '560px',
            h2("Optimal Number of Topics"),
            withSpinner(plotOutput("topics_plot"))
          ),
          box(
            width = 6,
            h2("Top Terms in Optimal Topic"),
            numericInput("num_top_terms", "Number of Top Terms to Display:", value = 3, min = 1, max = 3, step = 1),
            withSpinner(plotOutput("top_terms_plot"))
          )
        )
      ),
      # Differential Abundance tab
      tabItem(
        tabName = "diff_abundance",
        fluidRow(
          box(
            width = 12,
            h2("Differential Abundance of Topics"),
            withSpinner(plotOutput("diff_abundance_plot"))
          )
        )
      ),
      # Genera Analysis tab
      tabItem(
        tabName = "genera_analysis",
        fluidRow(
          box(
            width = 12,
            h2("Differential Abundance of Genera"),
            withSpinner(plotOutput("genera_plot"))
          )
        )
      ),
      # Selected Genera tab
      tabItem(
        tabName = "selected_genera",
        fluidRow(
          box(
            width = 6,
            h2("Analysis of Selected Genera"),
            withSpinner(plotOutput("selected_genera_plot"))
          ),
          box(
            width = 6,
            h2("Relative Abundance of Selected Genera"),
            withSpinner(plotOutput("relative_abundance_plot"))
          )
        )
      ),
      # About tab
      tabItem(
        tabName = "about",
        fluidRow(
          box(
            width = 12,
            h2("About This Dashboard"),
            p("This dashboard was created to analyze microbiome data and provide insights into the differential abundance of bacterial genera in colorectal cancer (CRC) patients compared to healthy controls."),
            h3("Data Source"),
            p("The data is sourced from the GuptaA_2019 study available in the curatedMetagenomicData package."),
            h3("How to Use"),
            p("Navigate through the tabs to explore different aspects of the analysis. Hover over plots for tooltips and additional information.")
          )
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Data preparation
  ps <- reactive({
    withProgress(message = 'Loading data...', value = 0, {
      # Increment progress
      incProgress(0.1, detail = "Loading metadata")
      meta_df <- curatedMetagenomicData::sampleMetadata
      
      incProgress(0.2, detail = "Filtering data")
      mydata <- meta_df %>%
        dplyr::filter(study_name == "GuptaA_2019") %>%
        dplyr::select(sample_id, subject_id, body_site, disease, age, gender, BMI) %>%
        dplyr::mutate(case_status = ifelse(disease == "CRC", "CRC",
                                           ifelse(disease == "healthy", "Control", "Other"))) %>%
        dplyr::filter(case_status != "Other") %>%
        na.omit()
      
      keep_id <- mydata$sample_id
      
      incProgress(0.3, detail = "Loading sample expressions")
      crc_se_sp <- sampleMetadata %>%
        dplyr::filter(sample_id %in% keep_id) %>%
        curatedMetagenomicData::returnSamples("relative_abundance", counts = TRUE)
      
      incProgress(0.4, detail = "Processing OTU table")
      crc_sp_df <- as.data.frame(as.matrix(assay(crc_se_sp)))
      
      crc_sp_df <- crc_sp_df %>%
        rownames_to_column(var = "Species")
      
      otu_tab <- crc_sp_df %>%
        tidyr::separate(Species, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
        dplyr::select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>%
        dplyr::mutate(Species = gsub("s__", "", Species)) %>%
        tibble::column_to_rownames(var = "Species")
      
      tax_tab <- crc_sp_df %>%
        dplyr::select(Species) %>%
        tidyr::separate(Species, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
        dplyr::mutate(Kingdom = gsub("k__", "", Kingdom),
                      Phylum = gsub("p__", "", Phylum),
                      Class = gsub("c__", "", Class),
                      Order = gsub("o__", "", Order),
                      Family = gsub("f__", "", Family),
                      Genus = gsub("g__", "", Genus),
                      Species = gsub("s__", "", Species)) %>%
        dplyr::mutate(spec_row = Species) %>%
        tibble::column_to_rownames(var = "spec_row")
      
      incProgress(0.6, detail = "Creating phyloseq object")
      rownames(mydata) <- NULL
      mydata <- mydata %>%
        tibble::column_to_rownames(var = "sample_id")
      
      ps <- phyloseq(sample_data(mydata),
                     otu_table(otu_tab, taxa_are_rows = TRUE),
                     tax_table(as.matrix(tax_tab)))
      
      sample_data(ps)$read_count <- sample_sums(ps)
      ps <- subset_samples(ps, read_count > 10^7)
      
      incProgress(0.8, detail = "Pruning taxa")
      minTotRelAbun <- 1e-5
      x <- taxa_sums(ps)
      keepTaxa <- (x / sum(x)) > minTotRelAbun
      ps <- prune_taxa(keepTaxa, ps)
      
      ps_g <- tax_glom(ps, taxrank = "Genus")
      taxa_names(ps_g) <- tax_table(ps_g)[, 6]
      
      incProgress(1, detail = "Data loading complete")
      return(ps_g)
    })
  })
  
  # Display sample data table using DT
  output$sample_table <- renderDataTable({
    sample_data(ps()) %>%
      as.data.frame() %>%
      rownames_to_column(var = "Sample ID")
  },
  options = list(
    pageLength = 10,
    scrollX = T,
    searching = FALSE,   # Disable global search box
    paging = F,       # Enable pagination
    info = FALSE,         # Hide table information summary
    dom = 't'
  ))
  
  # Plot case status distribution
  output$case_status_plot <- renderPlot({
    case_status_counts <- table(sample_data(ps())$case_status)
    barplot(case_status_counts, col = "steelblue", main = "Case Status Distribution", ylab = "Count")
  })
  
  # Optimal number of topics
  topics_result <- reactive({
    withProgress(message = 'Finding optimal number of topics...', value = 0, {
      set.seed(42)
      incProgress(0.2, detail = "Preparing count matrix")
      count_matrix <- data.frame(t(data.frame(otu_table(ps()))))
      
      # Remove terms and documents with zero counts
      count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
      count_matrix <- count_matrix[, colSums(count_matrix) > 0]
      
      incProgress(0.5, detail = "Running topic number optimization")
      result <- FindTopicsNumber(
        count_matrix,
        topics = seq(from = 3, to = 50, by = 3),
        metrics = c("CaoJuan2009", "Arun2010"),
        method = "VEM",
        control = list(seed = 42),
        mc.cores = 1,
        verbose = FALSE
      )
      
      incProgress(1, detail = "Optimization complete")
      return(result)
    })
  })
  
  output$topics_plot <- renderPlot({
    result <- topics_result()
    # Custom plotting function
    my_plot <- function (values) {
      if ("LDA_model" %in% names(values)) {
        values <- values[!names(values) %in% c("LDA_model")]
      }
      columns <- base::subset(values, select = 2:ncol(values))
      values <- base::data.frame(values["topics"], base::apply(columns,
                                                               2, function(column) {
                                                                 scales::rescale(column, to = c(0, 1), from = range(column))
                                                               }))
      values <- reshape2::melt(values, id.vars = "topics", na.rm = TRUE)
      values$group <- values$variable %in% c("Griffiths2004", "Deveaud2014")
      values$group <- base::factor(values$group, levels = c(FALSE, TRUE), labels = c("Minimize", "Maximize"))
      p <- ggplot(values, aes_string(x = "topics", y = "value", group = "variable")) +
      geom_line() + 
      geom_point(aes_string(shape = "variable"), size = 3) + 
      guides(size = FALSE, shape = guide_legend(title = "Metrics:")) + 
      scale_x_continuous(breaks = values$topics) + 
      labs(
        title = "Evaluation of Topic Models: Number of Topics vs Metrics",
        x = "\nNumber of Topics",
        y = NULL
      ) + 
      facet_grid(group ~ .) + 
      theme_bw() %+replace% theme(panel.grid.major.y = element_blank(),
                                           panel.grid.minor.y = element_blank(),
                                           panel.grid.major.x = element_line(colour = "grey70"),
                                           panel.grid.minor.x = element_blank(),
                                           legend.key = element_blank(),
                                           strip.text.y = element_text(angle = 90),
                                           legend.position = "bottom")
      return(p)
    }
    
    column_names <- c("topics", "CaoJuan2009", "Arun2010")
    p1 <- my_plot(result %>% select(all_of(column_names)))
    print(p1)
  })
  
  # Topic modeling and plotting top terms
  lda_results <- reactive({
    withProgress(message = 'Training LDA model...', value = 0, {
      incProgress(0.1, detail = "Calculating optimal topic number")
      result <- topics_result()
      
      # Optimal number of topics
      result$CaoJuan2009_norm <- (result$CaoJuan2009 - min(result$CaoJuan2009)) / (max(result$CaoJuan2009) - min(result$CaoJuan2009))
      result$Arun2010_norm <- (result$Arun2010 - min(result$Arun2010)) / (max(result$Arun2010) - min(result$Arun2010))
      result$mean_metric <- rowMeans(result[, c("CaoJuan2009_norm", "Arun2010_norm")])
      optimal_topic <- result[which.min(result$mean_metric), ]
      
      incProgress(0.3, detail = "Preparing count matrix")
      set.seed(42)
      count_matrix <- data.frame(t(data.frame(otu_table(ps()))))
      count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
      count_matrix <- count_matrix[, colSums(count_matrix) > 0]
      
      incProgress(0.6, detail = "Fitting LDA model")
      lda_k <- LDA(count_matrix, k = optimal_topic$topics, method = "VEM", control = list(seed = 42))
      
      incProgress(0.8, detail = "Extracting results")
      b_df <- data.frame(tidy(lda_k, matrix = "beta"))
      g_df <- data.frame(tidy(lda_k, matrix = "gamma")) %>%
        arrange(document, topic)
      
      incProgress(1, detail = "LDA training complete")
      list(b_df = b_df, g_df = g_df, optimal_topic = optimal_topic, lda_k = lda_k)
    })
  })
  
  output$top_terms_plot <- renderPlot({
    b_df <- lda_results()$b_df
    optimal_topic <- lda_results()$optimal_topic
    top1_topic <- which.min(lda_results()$optimal_topic$mean_metric)
    num_terms <- input$num_top_terms
    
    b_plot <- b_df %>%
      dplyr::filter(topic == num_terms) %>%
      arrange(desc(beta)) %>%
      slice_head(n = 20) %>%
      arrange(desc(beta)) %>%
      mutate(term = factor(term, levels = rev(term)))
    
    p3 <- ggplot(data = b_plot, aes(x = beta, y = term, color = beta)) +
      geom_point(aes(size = beta), alpha = 0.7) +
      scale_color_viridis_c(direction = -1) +
      labs(
        title = glue("Top 20 Terms in Topic {num_terms}"),
        subtitle = "Terms are ordered by their contribution to the topic",
        x = "\nBeta (Topic-Term Probability)",
        y = NULL,
        size = "Beta"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    print(p3)
  })
  
  # Differential abundance of topics
  output$diff_abundance_plot <- renderPlot({
    withProgress(message = 'Analyzing differential abundance...', value = 0, {
      incProgress(0.2, detail = "Preparing data")
      ps_topic_g <- {
        lda_k <- lda_results()$lda_k
        g_df <- lda_results()$g_df
        optimal_topic <- lda_results()$optimal_topic
        
        # Corrected line: Directly naming the column as 'read_count'
        lib_size_df <- data.frame(read_count = sample_sums(ps())) %>%
          rownames_to_column(var = "document")
        
        tm_df <- left_join(lib_size_df, g_df) %>%
          mutate(topic_count = read_count * gamma,
                 topic_count = round(topic_count, 0)) %>%
          dplyr::select(-read_count, -gamma) %>%
          pivot_wider(names_from = topic, values_from = topic_count) %>%
          dplyr::rename_with(~ paste0("Topic_", .), -document) %>%
          column_to_rownames(var = "document") %>%
          t(.) %>%
          data.frame(.)
        
        ps_topic_g <- phyloseq(
          sample_data(ps()),
          otu_table(tm_df, taxa_are_rows = TRUE))
      }
      
      incProgress(0.5, detail = "Running LinDA analysis")
      linda <- linda(otu.tab = data.frame(otu_table(ps_topic_g)), meta = data.frame(sample_data(ps_topic_g)),
                     formula = '~ case_status + age + gender + BMI', alpha = 0.05, winsor.quan = 0.97,
                     prev.cut = 0.3, lib.cut = 0, n.cores = 1)
      
      # Optionally, check the components of linda$output
      # print(names(linda$output))
      
      incProgress(0.8, detail = "Processing results")
      # Access the correct component
      linda_res_df <- data.frame(linda$output$case_status) %>%
        dplyr::arrange(padj) %>%
        rownames_to_column(var = "Topic")
      
      fdr_linda <- linda_res_df %>%
        mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
        separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
        mutate(Topic = gsub("_", " ", Topic)) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        mutate(Topic = factor(Topic, levels = rev(Topic)))
      
      incProgress(1, detail = "Plotting results")
      p2 <- ggplot(data = fdr_linda, aes(x = Topic, y = log2FoldChange, fill = reject)) +
        geom_col(alpha = .8) +
        labs(
          title = "Differential Abundance of Topics in CRC",
          y = "\nLog2 Fold-Change",
          x = ""
        ) +
        theme_minimal() + 
        theme(text= element_text(color = "black", size = 12,face = "bold"),
          legend.position = "none"
        ) +
        coord_flip() +
        scale_fill_manual(values = c("grey", "dodgerblue")) +
        geom_hline(yintercept = 0, linetype = "dotted")
      print(p2)
    })
  })
  
  
  
  # Differential abundance of genera
  output$genera_plot <- renderPlot({
    withProgress(message = 'Analyzing genera abundance...', value = 0, {
      incProgress(0.3, detail = "Running LinDA analysis on genera")
      linda_g_all <- linda(otu.tab = data.frame(otu_table(ps())), meta = data.frame(sample_data(ps())),
                           formula = '~ case_status + age + gender + BMI', alpha = 0.05, winsor.quan = 0.97,
                           prev.cut = 0, lib.cut = 0, n.cores = 1)
      
      linda_res_g_df <- data.frame(linda_g_all$output$case_statusCRC) %>%
        dplyr::arrange(padj) %>%
        rownames_to_column(var = "Genus")
      
      b_df <- lda_results()$b_df
      optimal_topic <- lda_results()$optimal_topic
      top1_topic <- which.min(lda_results()$optimal_topic$mean_metric)
      b_plot <- b_df %>%
        dplyr::filter(topic == optimal_topic$topics) %>%
        arrange(desc(beta)) %>%
        slice_head(n = 20) %>%
        arrange(desc(beta)) %>%
        mutate(term = factor(term, levels = rev(term)))
      keep_g <- b_plot$term
      
      incProgress(0.7, detail = "Processing results")
      log2_df <- linda_res_g_df %>%
        filter(Genus %in% keep_g) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        mutate(Genus = factor(Genus, levels = rev(Genus)))
      
      incProgress(1, detail = "Plotting results")
      p4 <- ggplot(data = log2_df, aes(x = Genus, y = log2FoldChange, fill = abs(log2FoldChange))) +
        geom_col(alpha = 0.8) +
        scale_fill_viridis_c(direction = -1) +
        labs(
          title = "Log2 Fold-Change in Genus Abundance",
          subtitle = "Comparing CRC vs Control",
          y = "Log2 Fold-Change",
          x = "",
          fill = "Log2 Fold-Change"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          legend.position = "none",
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5)
        ) +
        coord_flip() +
        geom_hline(yintercept = 0, linetype = "dotted", color = "red")
      print(p4)
    })
  })
  
  # Analysis of selected genera
  output$selected_genera_plot <- renderPlot({
    withProgress(message = 'Analyzing selected genera...', value = 0, {
      selected_genera <- c("Veillonella", "Flavonifractor", "Parabacteroides", "Alistipes", "Bifidobacterium")
      incProgress(0.3, detail = "Running LinDA analysis")
      linda_g_all <- linda(otu.tab = data.frame(otu_table(ps())), meta = data.frame(sample_data(ps())),
                           formula = '~ case_status + age + gender + BMI', alpha = 0.05, winsor.quan = 0.97,
                           prev.cut = 0, lib.cut = 0, n.cores = 1)
      
      linda_res_g_df <- data.frame(linda_g_all$output$case_statusCRC) %>%
        dplyr::arrange(padj) %>%
        rownames_to_column(var = "Genus")
      
      incProgress(0.6, detail = "Processing results")
      selected_genera_df <- linda_res_g_df %>%
        filter(Genus %in% selected_genera) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        mutate(Genus = factor(Genus, levels = rev(Genus)))
      
      incProgress(1, detail = "Plotting results")
      ggplot(data = selected_genera_df, aes(x = Genus, y = log2FoldChange, fill = log2FoldChange)) +
        geom_col(alpha = 0.8) +
        scale_fill_viridis_c(direction = -1) +
        labs(
          title = "Log2 Fold-Change in Selected Genus Abundance",
          subtitle = "Comparing CRC vs Control",
          y = "Log2 Fold-Change",
          x = "Bacterial Genus",
          fill = "Log2 Fold-Change"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(color = "black", size = 10),
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        ) +
        coord_flip()
    })
  })
  
  # Relative abundance of selected genera
  output$relative_abundance_plot <- renderPlot({
    selected_genera <- c("Veillonella", "Flavonifractor", "Parabacteroides", "Alistipes", "Bifidobacterium")
    otu_table_df <- data.frame(otu_table(ps()))
    sample_data_df <- data.frame(sample_data(ps())) %>%
      rownames_to_column(var = "Sample")
    
    selected_otu_table <- otu_table_df %>%
      rownames_to_column(var = "Genus") %>%
      filter(Genus %in% selected_genera) %>%
      gather(key = "Sample", value = "Abundance", -Genus)
    
    merged_data <- left_join(selected_otu_table, sample_data_df, by = "Sample")
    
    ggplot(data = merged_data, aes(x = Genus, y = Abundance, fill = case_status)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7, position = position_dodge(width = 0.75)) +
      scale_y_log10() +
      scale_fill_manual(values = c("dodgerblue", "coral"), labels = c("Control", "CRC")) +
      labs(
        title = "Relative Abundance of Selected Genera (Log Scale)",
        subtitle = "CRC vs Control Comparison",
        y = "Relative Abundance (Log10 Scale)",
        x = "Bacterial Genus",
        fill = "Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "top"
      ) +
      coord_flip()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
