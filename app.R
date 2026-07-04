library(shiny)
library(readr)
library(dplyr)
library(stringr)
library(visNetwork)
library(igraph)
library(bslib)

# Load data: prefer existing CSV in common relative locations; only fetch if missing
data_candidates <- c(
  file.path("data", "omnipath_interactions_table.csv"),
  file.path("..", "data", "omnipath_interactions_table.csv"),
  path.expand("~/Desktop/Science/Elena_R_scripts/retrieve-omnipath-app-main/data/omnipath_interactions_table.csv")
)
found <- data_candidates[file.exists(path.expand(data_candidates))]
if (length(found) > 0) {
  df <- readr::read_csv(found[1], show_col_types = FALSE)
} else if (exists("fetch_omnipath_data", mode = "function")) {
  df <- tryCatch(fetch_omnipath_data(cache_dir = "omnipathr-cache", out_file = "data/omnipath_interactions_table.csv"),
                 error = function(e) {
                   stop("Failed to fetch Omnipath data: ", conditionMessage(e))
                 })
  if (is.null(df) || (is.character(df) && !file.exists(df))) {
    stop("Omnipath data not available and fetch failed or returned nothing. Place data/omnipath_interactions_table.csv in the project.")
  }
} else {
  stop("Could not find data/omnipath_interactions_table.csv. Create the file or provide fetch_omnipath_data(). Checked: ", paste(data_candidates, collapse = ", "))
}

# basic checks and derive columns if needed
required_basic <- c("source_genesymbol", "target_genesymbol")
if (!all(required_basic %in% colnames(df))) {
  stop("CSV missing required columns: ", paste(setdiff(required_basic, colnames(df)), collapse = ", "))
}

if (!("interaction_type" %in% colnames(df))) {
  df <- df %>%
    dplyr::mutate(
      interaction_type = dplyr::case_when(
        (!is.na(consensus_stimulation) & as.numeric(consensus_stimulation) == 1) ~ "activation",
        (!is.na(consensus_inhibition) & as.numeric(consensus_inhibition) == 1) ~ "inhibition",
        TRUE ~ NA_character_
      )
    )
}
if (!("interaction_strength" %in% colnames(df))) {
  df <- df %>%
    dplyr::mutate(
      references = ifelse(is.na(references), "", as.character(references)),
      interaction_strength = stringr::str_split(references, pattern = "[;,|]") %>%
        lapply(function(x) sum(nzchar(trimws(x)))) %>%
        unlist() %>%
        as.numeric(),
      interaction_strength = ifelse(is.na(interaction_strength), 0, interaction_strength)
    )
}

df <- df %>%
  dplyr::mutate(
    interaction_strength = as.numeric(interaction_strength),
    interaction_type = as.character(interaction_type)
  )

genes <- sort(unique(na.omit(c(df$source_genesymbol, df$target_genesymbol))))

# UI
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly", base_font = font_google("Source Sans Pro")),
  tags$head(tags$style(HTML("
    :root {
      --pastel-primary: #A8D5E2;       /* pastel blue */
      --pastel-primary-hover: #86C9DB;
      --pastel-accent: #F7C6C7;        /* pastel coral for reset */
      --pastel-input: #F7FAFD;         /* pale input background */
      --pastel-text: #1b2430;
      --muted-border: rgba(27,36,48,0.08);
    }

    /* network labels and error color kept */
    .vis-network .node .label { font-weight: 600; }
    .shiny-output-error { color: #d9534f; }

    /* Buttons: primary (build) and generic buttons */
    .btn-primary, .btn {
      background-color: var(--pastel-primary) !important;
      border-color: var(--muted-border) !important;
      color: var(--pastel-text) !important;
      box-shadow: none !important;
      font-weight: 600;
    }
    .btn-primary:hover, .btn:hover {
      background-color: var(--pastel-primary-hover) !important;
      color: var(--pastel-text) !important;
    }

    /* View button (by id) styled with an accent pastel */
    #reset_view {
      background-color: var(--pastel-accent) !important;
      border-color: var(--muted-border) !important;
      color: var(--pastel-text) !important;
    }
    #reset_view:hover {
      background-color: #f4b3b4 !important;
    }

    /* Sidebar / panel background softened */
    .well, .sidebar .well {
      background: linear-gradient(180deg, #ffffff, #fbfeff) !important;
      border: 1px solid rgba(27,36,48,0.04) !important;
      box-shadow: none !important;
      border-radius: 8px;
    }

    /* Inputs (menus) pastel background */
    .selectize-control .selectize-input,
    .selectize-control .selectize-input input {
      background: var(--pastel-input) !important;
      border: 1px solid var(--muted-border) !important;
      color: var(--pastel-text) !important;
      border-radius: 6px !important;
    }
    .selectize-control.multi .selectize-input > div {
      background: transparent !important;
    }

    /* Contrast of selections */
    .selectize-control .option.active { background: rgba(134,201,219,0.4) !important; }
    .selectize-control .option { color: var(--pastel-text) !important; }

  "))),
  titlePanel(title = div(tags$strong("Retreive Your Knowledge Graph"))),
  sidebarLayout(
    sidebarPanel(
      helpText("Select source/target genes. Leave both empty to show top interactions."),
      selectizeInput("source_gene", "Source gene (multiple)", choices = NULL, multiple = TRUE,
                     options = list(placeholder = 'Type to search...', maxOptions = 1000)),
      selectizeInput("target_gene", "Target gene (multiple)", choices = NULL, multiple = TRUE,
                     options = list(placeholder = 'Type to search...', maxOptions = 1000)),
      numericInput("max_edges", "Max edges when showing defaults", value = 200, min = 10, step = 10),
      actionButton("build", "Build network", class = "btn-primary"),
      width = 3
    ),
    mainPanel(
      div(style = "display:flex; align-items:center; justify-content:space-between;",
          verbatimTextOutput("info"),
          actionButton("reset_view", "Reset view")
      ),
      visNetworkOutput("network", height = "720px"),
      hr(),
      tableOutput("interactions"),
      tags$div(style = "margin-top:12px; color:#6c757d; font-size:0.9em; text-align:right;",
           "Author: Elena Eyre-Sanchez, PhD — Last updated: 2026"),
      width = 9
    )
  )
)

# internal renderer function
server_network_render <- function(eds) {
  if (nrow(eds) == 0) {
    return(visNetwork::visNetwork(data.frame(id = character(0)), data.frame()))
  }
  
  nodes_vec <- sort(unique(c(eds$source_genesymbol, eds$target_genesymbol)))
  nodes <- data.frame(id = nodes_vec, label = nodes_vec, title = nodes_vec, stringsAsFactors = FALSE)
  
  g <- igraph::graph_from_data_frame(dplyr::select(eds, source_genesymbol, target_genesymbol),
                                     directed = TRUE, vertices = nodes)
  deg <- igraph::degree(g, mode = "all")
  deg_vec <- setNames(as.numeric(deg), names(deg))
  nodes$size <- 12 + 4 * (ifelse(is.na(deg_vec[as.character(nodes$id)]), 0, deg_vec[as.character(nodes$id)]))
  nodes$shadow <- TRUE
  
  # map types to colors
  palette_fixed <- c(activation = "#b31b1c", inhibition = "#e31a1c", unknown = "#6a6a6a")
  types <- unique(na.omit(eds$interaction_type))
  if (length(types) == 0) types <- "unknown"
  type_colors <- palette_fixed[names(palette_fixed) %in% types]
  missing_types <- setdiff(types, names(type_colors))
  if (length(missing_types)) {
    extra_cols <- grDevices::colorRampPalette(c("#8dd3c7","#ffffb3","#bebada"))(length(missing_types))
    type_colors <- c(type_colors, setNames(extra_cols, missing_types))
  }
  edge_colors <- ifelse(is.na(eds$interaction_type), "#999999", type_colors[eds$interaction_type])
  
  max_abs <- suppressWarnings(max(abs(as.numeric(eds$interaction_strength)), na.rm = TRUE))
  if (is.na(max_abs) || max_abs == 0) max_abs <- 1

  # build edge dataframe for visNetwork (do not include a single-element list column)
  eds_vis <- data.frame(
    from = eds$source_genesymbol,
    to   = eds$target_genesymbol,
    value = pmax(1, 6 * (abs(as.numeric(eds$interaction_strength)) / max_abs)),
    title = paste0("<b>Type:</b> ", ifelse(is.na(eds$interaction_type),"NA",eds$interaction_type),
                   "<br><b>Strength:</b> ", eds$interaction_strength),
    color = edge_colors,
    arrows = "to",
    stringsAsFactors = FALSE
  )
  
  vis <- visNetwork(nodes, eds_vis) %>%
    visEdges(color = list(highlight = "#ffeb3b"), smooth = TRUE) %>%
    visNodes(shape = "dot", shadow = TRUE, font = list(size = 14, face = 'sans')) %>%
    visLayout(randomSeed = 123) %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
               nodesIdSelection = list(enabled = TRUE, useLabels = TRUE)) %>%
    visPhysics(stabilization = list(enabled = TRUE, iterations = 1000),
               barnesHut = list(gravitationalConstant = -30000, springLength = 200)) %>%
    visLegend(useGroups = FALSE,
              addNodes = lapply(names(type_colors), function(t) list(label = t, shape = "dot", color = type_colors[[t]], size = 12)),
              width = 0.2, position = "right")
  vis
}

# Server
server <- function(input, output, session) {
  # populate selectize with server = TRUE
  updateSelectizeInput(session, "source_gene", choices = genes, server = TRUE)
  updateSelectizeInput(session, "target_gene", choices = genes, server = TRUE)
  
  build_edges <- eventReactive(input$build, {
    req(df)
    res <- df
    
    if (!is.null(input$source_gene) && length(input$source_gene) > 0) {
      res <- dplyr::filter(res, source_genesymbol %in% input$source_gene)
    }
    if (!is.null(input$target_gene) && length(input$target_gene) > 0) {
      res <- dplyr::filter(res, target_genesymbol %in% input$target_gene)
    }
    
    if (nrow(res) == 0) {
      showNotification("No edges match — showing top interactions by strength.", type = "warning")
      res <- df %>%
        dplyr::mutate(abs_strength = abs(as.numeric(interaction_strength))) %>%
        dplyr::arrange(desc(abs_strength)) %>%
        dplyr::slice_head(n = input$max_edges) %>%
        dplyr::select(-abs_strength)
    } else {
      if (nrow(res) > 2000) {
        showNotification("Filtered result large — showing strongest 2000 edges.", type = "message")
        res <- res %>%
          dplyr::mutate(abs_strength = abs(as.numeric(interaction_strength))) %>%
          dplyr::arrange(desc(abs_strength)) %>%
          dplyr::slice_head(n = 2000) %>%
          dplyr::select(-abs_strength)
      }
    }
    
    res <- dplyr::mutate(res, interaction_strength = as.numeric(interaction_strength))
    res
  }, ignoreNULL = FALSE)
  
  output$interactions <- renderTable({
    head(build_edges(), 200)
  })
  
  output$info <- renderText({
    eds <- build_edges()
    paste0("Edges: ", nrow(eds), " | Unique genes: ", length(unique(c(eds$source_genesymbol, eds$target_genesymbol))))
  })
  
  # use the internal renderer here
  output$network <- renderVisNetwork({
    server_network_render(build_edges())
  })
  
  # reset view button: fit the network
  observeEvent(input$reset_view, {
    try({
      visNetworkProxy("network") %>% visFit()
    }, silent = TRUE)
  })
}

# sanity check and launch
if (!exists("server") || !is.function(server)) {
  stop("The Shiny 'server' function is not defined. Check app.R for syntax errors.")
}
if (!exists("ui")) {
  stop("The Shiny 'ui' object is not defined.")
}

shiny::shinyApp(ui = ui, server = server)
