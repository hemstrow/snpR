library(shiny)
source("app_code.R")

# other things I'd like this to do:
## citations
## monitor and save repeatable snpR code

ui <- fluidPage(
  
  #===============input===============
  titlePanel("Data import:"),
  sidebarLayout( 
    sidebarPanel(
      fluidRow(fileInput(inputId = "genotypes", label = "genotypes/vcf/ms/etc"),
               fileInput(inputId = "snp_meta", label = "SNP metadata (optional)"),
               fileInput(inputId = "sample_meta", label = "sample metadata (optional)"),
               textInput(inputId = "mDat", label = "missing data format (optional)", value = "NN", placeholder = "NN, 0000, etc."),
               actionButton("import", "Import Data"))
    ),
    mainPanel(fluidRow(
      plotOutput("input_validation"),
      textOutput("input_validation_numbers")))
  ),
  
  # subset sidebar
  textOutput(outputId = "subset_header"),
  fluidRow(textInput(inputId = "subset", label = "Subsetting call, see ?subset_snpR_data for more info"),
           actionButton("subset", "Subset Data")), # note that this should change the input_validation return!

  #============filtering================
  titlePanel("Filtering:"),
  fluidRow(
    column(5,
           titlePanel(h4("Filtering option A:")),
           numericInput(inputId = "mafA", label = "Minor Allele Frequency Minimum", value = 0, min = 0, max = .5),
           selectInput("maf_facets_A", "MAF facet: keep any loci above specified Minor Alele Frequency in any facet level.", choices = list("Upload data to see available facets." = "none"),
                       multiple = TRUE),
           numericInput(inputId = "hweA", label = "HWE p-value Maximum", value = 0, min = 0, max = 1),
           selectInput("hwe_facets_A", "HWE facet: reject any loci out of HWE in any facet level", choices = list("Upload data to see available facets." = "none"),
                       multiple = TRUE),
           numericInput(inputId = "min_indA", label = "Maximum Proportion of Missing SNPs per Sample", value = 0, min = 0, max = 1),
           numericInput(inputId = "min_lociA", label = "Maximum Proportion of Missing Samples per Locus", value = 0, min = 0, max = 1),
           actionButton("run_filt_A", "Run?")),
    column(5,
           titlePanel(h4("Filtering option B:")),
           numericInput(inputId = "mafB", label = "Minor Allele Frequency Minimum", value = 0),
           selectInput("maf_facets_B", "MAF facet: keep any loci above specified Minor Alele Frequency in any facet level.", choices = list("Upload data to see available facets." = "none"),
                       multiple = TRUE),
           numericInput(inputId = "hweB", label = "HWE p-value Maximum", value = 0),
           selectInput("hwe_facets_B", "HWE facet: reject any loci out of HWE in any facet level", choices = list("Upload data to see available facets." = "none"),
                       multiple = TRUE),
           numericInput(inputId = "min_indB", label = "Maximum Proportion of Missing SNPs per Sample", value = 0),
           numericInput(inputId = "min_lociB", label = "Maximum Proportion of Missing Samples per Locus", value = 0),
           actionButton("run_filt_B", "Run?"))
  ),
  fluidRow(column(5,
                  fluidRow(
                    plotOutput("filter_validation_A"),
                    textOutput("filter_validation_A_numbers")
                  ),
                  actionButton("filter_select_A", "Apply filtering option A")),
           column(5,
                  fluidRow(
                    plotOutput("filter_validation_B"),
                    textOutput("filter_validation_B_numbers")
                  ),
                  actionButton("filter_select_B", "Apply filtering option B"))
  ),

  #===========statistics=======



  # tests - single
  uiOutput(outputId = "stats_header"),
  checkboxGroupInput("selected_tests_single", label = "Desired Single-SNP Statistics:",
                     choices = c("Ho" = "ho",
                                 "pi" = "pi",
                                 "Fst" = "Fst",
                                 "Fis" = "Fis",
                                 "HWE" = "hwe",
                                 "Minor Allele Frequencies" = "maf",
                                 "Private Alleles" = "pa")),
  fluidRow(checkboxInput("do_fst_boot", "Calculate Fst p-values via Bootstrapping?")),
  textInput("fst_boots", "# Bootstrapps", value = 1000),
  uiOutput("single_facets"), # provide drop down with clickable facets to use

  ## tests - window
  checkboxGroupInput("selected_tests_window", label = "Desired Windowed Statistics",
                     choices = c("Ho" = "ho",
                                 "pi" = "pi",
                                 "Fst" = "Fst",
                                 "Private Alleles" = "pa",
                                 "Tajima's D" = "tsd")),
  uiOutput(outputId = "window_opts"),
  fluidRow(textInput(inputId = "sigma", label = "Sigma (Window Size = 6*Sigma, in kb)", value = 200),
           textInput(inputId = "slide", label = "Slide Between Windows (kb)", value = 50),
           checkboxInput(inputId = "do_boots", "Bootstrap Window Significance?"),
           textInput(inputId = "num_window_boots", "Number of Bootstraps", value = 1e6)),
  uiOutput("window_facets"),

  ## association

  ## plotting

  ## tests - misc
  fluidRow(
    uiOutput("sfs_facets"),
    uiOutput("sfs_facet_levels")
  ),
  fluidRow(checkboxInput("ibd", "Isolation by Distance"),
           uiOutput("snp_ibd_facets")),
  fluidRow(checkboxInput("ne", "Effective Population Size (Ne)"),
           textInput("neestimator_path", "Path to NeEstimator executable"),
           checkboxGroupInput("ne_methods", label = "Methods:", choices = list("LDNe" = "ld",
                                                                               "Heterozygote Excess" = "Ht",
                                                                               "Coancestry" = "coan")),
           numericInput(inputId = "ne_pcrit", label = "P-crit (minimum minor allele frequency)", 0.01, 0, .5))
)



server <- function(input, output, session) {
  
  #==========read in input files===========
  x <- reactiveValues()

  data_importer <- reactive({
    x <- import_wrapper(input$genotypes, input$snp_meta, input$sample_meta, input$mDat)
    valid <- validation(x)
    return(list(x = x, valid = valid))
  })
  
  observeEvent(input$import,{
    x$data <- data_importer()
  })
  
  
  observeEvent(input$import | input$filter_select_A | input$filter_select_B, ignoreInit = TRUE, {
    cat("updated plot")
    vpns <- renderValidation(x$data$valid)
    output$input_validation <- vpns$vp
    output$input_validation_numbers <- vpns$vn
  })


  #==========filtering=====================
  # set up reactives
  filtered_sets <- reactiveValues()
  a <- reactive({
    d <- filter_wrapper(x$data$x,
                        maf = input$mafA,
                        maf_facet = input$maf_facets_A,
                        hwe = input$hweA,
                        hwe_facet = input$hwe_facets_A,
                        min_ind = input$min_indA,
                        min_loci = input$min_lociA)
    
    v <- validation(d)
    return(list(x = d, valid = v))
  })
  b <- reactive({
    d <- filter_wrapper(x$data$x,
                        maf = input$mafB,
                        maf_facet = input$maf_facets_B,
                        hwe = input$hweB,
                        hwe_facet = input$hwe_facets_B,
                        min_ind = input$min_indB,
                        min_loci = input$min_indB)
    
    v <- validation(d)
    return(list(x = d, valid = v))
  })
  
  # facet level trackers
  observe({
    
    if(!is.null(x$data$x)){
      facet_opts <- facet_check(x$data$x)
      
      updateSelectInput(session, "maf_facets_A", choices = facet_opts$sample)
      updateSelectInput(session, "maf_facets_B", choices = facet_opts$sample)
      updateSelectInput(session, "hwe_facets_A", choices = facet_opts$sample)
      updateSelectInput(session, "hwe_facets_B", choices = facet_opts$sample)
    }
    
    
  })
  
  # plot filters
  observeEvent(input$run_filt_A,{
    filtered_sets$a <- a()
    av <- validation(filtered_sets$a$x)
      
    vpnsa <- renderValidation(av)
    output$filter_validation_A <- vpnsa$vp
    output$filter_validation_A_numbers <- vpnsa$vn
  })
  
  observeEvent(input$run_filt_B,{
    filtered_sets$b <- b()
    bv <- validation(filtered_sets$b$x)
    
    vpnsb <- renderValidation(bv)
    output$filter_validation_B <- vpnsb$vp
    output$filter_validation_B_numbers <- vpnsb$vn
  })  
  
  # keep a filter
  observeEvent(input$filter_select_A, {
    x$data <- filtered_sets$a
  })
  observeEvent(input$filter_select_B, {
    x$data <- filtered_sets$b
  })

  #===================do statistics=====================
  
} 

shinyApp(ui, server)