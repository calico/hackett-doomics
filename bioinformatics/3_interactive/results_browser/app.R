library(dplyr)
library(ggplot2)
library(patchwork)
library(shiny)

source("interactive.R")

# read results from GCP
domics_paths <- read_gcp_modeling_results(outdir = "/tmp", overwrite = FALSE)
model_signif <- readRDS(domics_paths$local_path[domics_paths$type == "signif"])
features_with_design <- readRDS(domics_paths$local_path[domics_paths$type == "features"])

# populate data modalities
data_modalities <- unique(model_signif$data_modality)
all_features <- model_signif %>%
  dplyr::distinct(data_modality, groupName)

ui = fluidPage(
  tags$head(tags$style(
    type = "text/css",
    "h1, h2, h3, h4, h5, h6 { color: #5BB867;}",
    "label { font-size: 20px;}",
    "div { font-size: 15px;}",
    "body {width: 100% !important; max-width: 100% !important;}"
  )),
  
  # Application title
  headerPanel("DO Results Browser"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      shiny::h3("Select Data Modalities"),
      wellPanel(
        selectizeInput(
          "selected_modalities",
          "Data Modalities:",
          choices = data_modalities,
          selected = data_modalities,
          multiple = TRUE
        )
      ),
      shiny::h3("Select Molecules"),
      wellPanel(
        selectizeInput(
          "selected_features",
          "Molecules:",
          choices = NULL,
          multiple = TRUE
        )
      )
    ),
    mainPanel(
      shiny::tabsetPanel(
        type = "tabs",
        shiny::tabPanel(
          "Summary",
          shiny::dataTableOutput("overall_signif")
        ),
        shiny::tabPanel("Sex", signifOutput("sex")),
        shiny::tabPanel("Age", signifOutput("aging")),
        shiny::tabPanel("DDM", signifOutput("ddm")),
        shiny::tabPanel("Lifespan", signifOutput("lifespan")),
        shiny::tabPanel("Lifespan Remaining", signifOutput("lifespan_remaining")),
        shiny::tabPanel("Fraction of Life Lived", signifOutput("fraction_of_life_lived")),
        shiny::tabPanel("Age x Lifespan", signifOutput("age_by_lifespan"))
      )
    )
  )
)

server = function(input, output, session) {
  
  # reduce data and features based on selected data modalities
  
  selected_datasets <- reactive({
    req(input$selected_modalities)
    
    features_with_design %>%
      dplyr::filter(data_modality %in% input$selected_modalities)
  })
  
  selected_signif <- reactive({
    req(input$selected_modalities)
    
    model_signif %>%
      dplyr::filter(data_modality %in% input$selected_modalities)
  })
  
  observe({
    req(selected_signif())
    
    all_available_features <- all_features %>%
      dplyr::filter(data_modality %in% input$selected_modalities)
    
    updateSelectizeInput(
      session,
      "selected_features",
      choices = all_available_features$groupName,
      server = TRUE
    )
  })
  
  # reduce data based on selected features
  
  selected_features_signif <- reactive({
    req(input$selected_features)
    
    selected_signif() %>%
      dplyr::filter(groupName %in% input$selected_features)
  })
  
  
  selected_features_data <- reactive({
    req(input$selected_features)
    
    selected_datasets() %>%
      dplyr::semi_join(
        selected_features_signif(),
        by = c("data_type", "groupId")
      )
  })

  observe({
    
    # all tests of features
    output$overall_signif <- shiny::renderDataTable(
      selected_features_signif() %>%
        clean_plot_signif() %>%
        dplyr::arrange(q.value)
        )
    
    signifServer("aging", plot_aging(selected_features_data(), selected_features_signif(), input$selected_features))
    signifServer("lifespan", plot_lifespan(selected_features_data(), selected_features_signif(), input$selected_features))
    signifServer("lifespan_remaining", plot_lifespan_remaining(selected_features_data(), selected_features_signif(), input$selected_features))
    signifServer("age_by_lifespan", plot_age_by_lifespan(selected_features_data(), selected_features_signif(), input$selected_features))
    signifServer("ddm", plot_ddm(selected_features_data(), selected_features_signif(), input$selected_features))
    signifServer("sex", plot_sex(selected_features_data(), selected_features_signif(), input$selected_features))
    signifServer("fraction_of_life_lived", plot_fraction_of_life_lived(selected_features_data(), selected_features_signif(), input$selected_features))
  })
}

shiny::shinyApp(ui = ui, server = server)
