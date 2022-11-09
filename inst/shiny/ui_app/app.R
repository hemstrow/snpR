library(shiny)
library(shinyjs)
source("app_code.R")

# other things I'd like this to do:
## citations
## monitor and save repeatable snpR code

ui <- navbarPage(
  
  #===============input===============
  tabPanel("Data Import",
           useShinyjs(),
           titlePanel("Data import:"),
           actionButton("use_example", "Use example data?"),
           sidebarLayout( 
             sidebarPanel(
               fluidRow(fileInput(inputId = "genotypes", label = "genotypes/vcf/ms/etc"),
                        fileInput(inputId = "snp_meta", label = "SNP metadata (optional)"),
                        fileInput(inputId = "sample_meta", label = "sample metadata (optional)"),
                        textInput(inputId = "mDat", label = "missing data format (optional)", value = "NN", placeholder = "NN, 0000, etc."),
                        uiOutput("import"))
             ),
             mainPanel(fluidRow(
               plotOutput("input_validation"),
               textOutput("input_validation_numbers")))
           ),
           
           # subset sidebar
           textOutput(outputId = "subset_header"),
           fluidRow(textInput("subsetting_call", "Subsetting call, see ?subset_snpR_data for more info"),
                    actionButton("subset", "Subset Data")) # note that this should change the input_validation return!
  ),
  #============filtering================

  tabPanel("Filtering",
           titlePanel("Filtering:"),
           fluidRow(
             column(6,
                    titlePanel(h4("Filtering option A:")),
                    numericInput(inputId = "mafA", label = "Minor Allele Frequency Minimum", value = 0, min = 0.05, max = .5),
                    selectInput("maf_facets_A", "MAF facet: keep any loci above specified Minor Alele Frequency in any facet level.", choices = list("Upload data to see available facets." = "none"),
                                multiple = TRUE),
                    checkboxInput("singletons_A", "Remove singletons"),
                    numericInput(inputId = "hweA", label = "HWE p-value Maximum", value = 0.0001, min = 0, max = 1),
                    selectInput("hwe_facets_A", "HWE facet: reject any loci out of HWE in any facet level", choices = list("Upload data to see available facets." = "none"),
                                multiple = TRUE),
                    numericInput(inputId = "min_indA", label = "Maximum Proportion of Missing SNPs per Sample", value = 0.75, min = 0, max = 1),
                    numericInput(inputId = "min_lociA", label = "Maximum Proportion of Missing Samples per Locus", value = 0.75, min = 0, max = 1),
                    actionButton("run_filt_A", "Run?")),
             column(6,
                    titlePanel(h4("Filtering option B:")),
                    numericInput(inputId = "mafB", label = "Minor Allele Frequency Minimum", value = 0),
                    selectInput("maf_facets_B", "MAF facet: keep any loci above specified Minor Alele Frequency in any facet level.", choices = list("Upload data to see available facets." = "none"),
                                multiple = TRUE),
                    checkboxInput("singletons_B", "Remove singletons"),
                    numericInput(inputId = "hweB", label = "HWE p-value Maximum", value = 0.000001),
                    selectInput("hwe_facets_B", "HWE facet: reject any loci out of HWE in any facet level", choices = list("Upload data to see available facets." = "none"),
                                multiple = TRUE),
                    numericInput(inputId = "min_indB", label = "Maximum Proportion of Missing SNPs per Sample", value = 0.5, min = 0, max = 1),
                    numericInput(inputId = "min_lociB", label = "Maximum Proportion of Missing Samples per Locus", value = 0.5, min = 0, max = 1),
                    actionButton("run_filt_B", "Run?"))
           ),
           fluidRow(column(6,
                           fluidRow(
                             plotOutput("filter_validation_A"),
                             textOutput("filter_validation_A_numbers")
                           ),
                           actionButton("filter_select_A", "Apply filtering option A")),
                    column(6,
                           fluidRow(
                             plotOutput("filter_validation_B"),
                             textOutput("filter_validation_B_numbers")
                           ),
                           actionButton("filter_select_B", "Apply filtering option B"))
           )
  ),

  #===========statistics=======
  tabPanel("Statistics",
           tabsetPanel(


             # tests - single
             tabPanel("Basic Analysis",
                      uiOutput(outputId = "stats_header"),
                      checkboxGroupInput("selected_tests_single", label = "Options:",
                                         choices = c("Ho" = "ho",
                                                     "He" = "he",
                                                     "pi" = "pi",
                                                     "Fst" = "Fst",
                                                     "Fis" = "Fis",
                                                     "HWE" = "hwe",
                                                     "Hs (individual heterozygosity)" = "Hs",
                                                     "Minor Allele Frequencies" = "maf",
                                                     "Private Alleles" = "pa")),
                      checkboxInput("do_fst_boot", "Calculate Fst p-values via Bootstrapping?"),
                      hidden(textInput("fst_boots", "# Bootstrapps", value = 1000)),
                      uiOutput("single_facets"), # provide drop down with clickable facets to use
                      actionButton("run_single_tests", "Run!")
             ),


             # tests - window
             tabPanel("Sliding Window Analysis",
                      checkboxGroupInput("selected_tests_window", label = "Options:",
                                         choices = c("Ho" = "ho",
                                                     "pi" = "pi",
                                                     "Fst" = "Fst",
                                                     "Private Alleles" = "pa",
                                                     "Tajima's D" = "tsd")),
                      uiOutput(outputId = "window_opts"),
                      fluidRow(textInput(inputId = "sigma", label = "Sigma (Window Size = 6*Sigma, in kb)", value = 200),
                               textInput(inputId = "slide", label = "Slide Between Windows (kb)", value = 50),
                               checkboxInput(inputId = "do_boots", "Bootstrap Window Significance?"),
                               hidden(textInput(inputId = "num_window_boots", "Number of Bootstraps", value = 1e6))),
                      uiOutput("window_facets"),
                      numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1),
                      actionButton("run_window_tests", "Run!")
             ),

             # tests - association
             tabPanel("Association Testing and Genomic Prediction",
                      radioButtons("association_method", "Method:", choices = list(GWAS = "GWAS",
                                                                                   GP = "Genomic Prediction",
                                                                                   RF = "Random Forest"), inline = TRUE),
                      uiOutput("all_facets"),
                      uiOutput("association_response"),
                      uiOutput("association_covariates"),
                      hidden(tags$div(id = "GWAS",
                                      selectInput("GWAS_method", label = "Method:", choices = c("gmmat.score" = "GMMAT",
                                                                                                "armitage" = "Armitage",
                                                                                                "odds_ratio" = "Log Odds Ratio",
                                                                                                "chisq" = "Chi-squared")),
                                      hidden(tags$div(id = "armitage_options",
                                                      inputPanel(
                                                        fluidRow(h5("Armitage Test Weights:"),
                                                                 numericInput("w1", "Homozygote A", 0, 0, 1),
                                                                 numericInput("w2", "Heterozygote", .5, 0, 1),
                                                                 numericInput("w2", "Homozygote B", 1, 0, 1))))),
                                      hidden(tags$div(id = "GMMAT_options",
                                                      inputPanel(
                                                        uiOutput("GMMAT_covariates"),
                                                        textInput("GWAS_family", "Regression family override, must be a valid family. See ?family."),
                                                        numericInput("maxiter", "Maximum number of fitting iterations", 500, 100, 10000),
                                                        uiOutput("GMMAT_sampleID"),
                                                        numericInput("Gmaf", "Minimum Minor Allele Frequency?", 0, 0, .5),
                                                        numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1)))))),

                      hidden(tags$div(id = "GP",
                                      numericInput("GP_iters", 10000, 100, 1e7),
                                      uiOutput("GP_burnin"),
                                      uiOutput("GP_thin"),
                                      selectInput("GP_model", "Genetic Architecture Model:",
                                                  choices = c("BayesA", "BayesB", "BayesC", "BRR", "BL", "FIXED"), selected = "BayesB"),
                                      selectInput("GP_interpolate", "Missing Genotype Interpolation",
                                                  choices = c("binomial draw" = "bernoulli", "allele frequency" = "af", "iPCA (warning: slow)" = "iPCA")),
                                      fluidRow(hidden(numericInput("GP_iPCA_ncp", "ncp:", NULL, 0)), hidden(numericInput("GP_iPCA_ncp_max", "ncp max:", 1, 5))),
                                      numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1))),


                      hidden(tags$div(id = "RF",
                                      uiOutput("RF_covariates"),
                                      numericInput("RF_num.trees", "Number of Trees", 10000, 100, 1e7),
                                      uiOutput("RF_mtry"),
                                      selectInput("RF_importance", "Importance method:", choices = c("", impurity_corrected = "Corrected Impurity", impurity = "Impurity", permutation = "Permutation")),
                                      selectInput("RF_interpolate", "Missing Genotype Interpolation",
                                                  choices = c(`binomial draw` = "bernoulli", `allele frequency` = "af", `iPCA (warning: slow)` = "iPCA")),
                                      fluidRow(hidden(numericInput("iPCA_ncp", "ncp:", NULL, 0)), hidden(numericInput("iPCA_ncp_max", "ncp max:", 1, 5))),
                                      numericInput("par", label = "Number of parallel processing threads:", value = 1, min = 1, max = parallel::detectCores(), step = 1)))


             ),

             ## tests - misc
             tabPanel("Misc",
                      fluidRow(
                        uiOutput("sfs_facets"),
                        uiOutput("sfs_facet_levels")
                      ),
                      fluidRow(checkboxInput("ibd", "Isolation by Distance"),
                               uiOutput("snp_ibd_facets")),
                      fluidRow(checkboxInput("ne", "Effective Population Size (Ne)"),
                               textInput("neestimator_path", "Path to NeEstimator executable", "/usr/bin/Ne2-1.exe"),
                               checkboxGroupInput("ne_methods", label = "Methods:", choices = list("LDNe" = "ld",
                                                                                                   "Heterozygote Excess" = "Ht",
                                                                                                   "Coancestry" = "coan")),
                               numericInput(inputId = "ne_pcrit", label = "P-crit (minimum minor allele frequency)", 0.01, 0, .5))
             )


           )
  ),

  #================plotting=================
  # tabPanel("Plotting",
  #          # PCA
  #          tabsetPanel(
  #            tabPanel("PCA/UMAP/tSNE",
  #                     uiOutput("all_facets"),
  #                     radioButtons("association_method", "Method:", choices = list(PCA = "pca",
  #                                                                                  UMAP = "umap",
  #                                                                                  tSNE = "tsne"), inline = TRUE),
  #                     checkboxInput("check_duplicates", "Check for Duplicates (slow)?", FALSE),
  #                     numericInput("minimum_percent_coverage", "Minimum proportion of genotypes per sample:", 0, 0, 1),
  #                     numericInput("minimum_genotype_percentage", "Minimum proportion of samples sequenced per loci:", 0, 0, 1),
  #                     selectInput("clusters_interpolate", "Missing Genotype Interpolation",
  #                                 choices = c(bernoulli = "binomial draw", af = "allele frequency", iPCA = "iPCA (warning: slow)")),
  #                     fluidRow(hidden(numericInput("iPCA_ncp", "ncp:", NA, 0)), hidden(numericInput("iPCA_ncp_max", "ncp max:", 5, 1))),
  #                     hidden(tags$div(id = "tSNE_options",
  #                                     inputPanel(numericInput("tSNE_dims", "Output dimensions (at max 2 plotted):", 2, 1),
  #                                                numericInput("tSNE_perplexity", "Perplexity:", NA, min = 0))))
  #            ),
  # 
  #            # structure
  #            tabPanel("STRUCTURE/ADMIXTURE q-plots",
  #                     uiOutput("all_facets"),
  #                     uiOutput("facet_order"),
  #                     numericInput("kmin", "Minimum k value:", 1, 1),
  #                     numericInput("kmin", "Maximum k value:", 2, 1),
  #                     selectInput("struc_method", "Assignment/Clustering Method:",
  #                                 choices = c(STRUCTURE = "structure",
  #                                             sNMF = "snmf",
  #                                             ADMIXTURE = "admixture",
  #                                             snapclust = "snapclust"),
  #                                 selected = "snmf")
  #            )
  #          )
  #          # LD
  #          # manhattan
  #          # sfs
  # ),

  #================navbar options===========
  title = "snpR"
  
  #================end of ui================
)



server <- function(input, output, session) {

  #==========read in input files===========
  x <- reactiveValues(x = NULL,
                      valid = NULL)

  data_importer <- reactive({
    return(import_wrapper(input$genotypes, input$snp_meta, input$sample_meta, input$mDat))
  })
  
  observeEvent(input$import,{
    x$x <- data_importer()
  })
  
  
  observe({
    if(is.snpRdata(x$x)){
      isolate(x$valid <- validation(x$x))
    }
  })
  
  
  observeEvent(input$use_example, ignoreInit = TRUE, {
    ex <- stickSNPs
    ex <- calc_fis(ex)
    x$x <- ex
  })
  
  
  observeEvent(x$valid, ignoreInit = TRUE, {
    cat("updated plot")
    if(is.list(x$valid)){
      
      vpns <- renderValidation(x$valid)
      
      output$input_validation <- vpns$vp
      output$input_validation_numbers <- vpns$vn
    }
  })
  
  # define import button
  output$import <- renderUI({
    validate(need(is.data.frame(input$genotypes), "Select a genotypes file to continue!"))
    
    actionButton("import", "Import Data")
  })
  
  # subsetting
  # observeEvent(input$subset,{
  #   x$x <- subset_wrapper(x$x, input$subsetting_call)
  # })


  #==========filtering=====================
  # set up reactives
  filtered_sets <- reactiveValues()
  a <- reactive({
    d <- filter_wrapper(x$x,
                        maf = input$mafA,
                        maf_facet = input$maf_facets_A,
                        singletons = input$singletons_A,
                        hwe = input$hweA,
                        hwe_facet = input$hwe_facets_A,
                        min_ind = input$min_indA,
                        min_loci = input$min_lociA)
    
    v <- validation(d)
    return(list(x = d, valid = v))
  })
  b <- reactive({
    d <- filter_wrapper(x$x,
                        maf = input$mafB,
                        maf_facet = input$maf_facets_B,
                        singletons = input$singletons_B,
                        hwe = input$hweB,
                        hwe_facet = input$hwe_facets_B,
                        min_ind = input$min_indB,
                        min_loci = input$min_indB)
    
    v <- validation(d)
    return(list(x = d, valid = v))
  })
  
  # facet level trackers
  observe({
    
    if(is.snpRdata(x$x)){
      facet_opts <- facet_check(x$x)
      
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
    x$x <- filtered_sets$a
  })
  observeEvent(input$filter_select_B, {
    x$x <- filtered_sets$b
  })

  #===================do statistics=====================
  
} 

shinyApp(ui, server)