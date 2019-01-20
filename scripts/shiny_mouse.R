library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)

WINDOW_SIZE = "700px"

tsne_dat = readRDS("DrugEffects_AllDrugs_EffectLength_tSNE_mouse_25components.rds")

# A bug in shiny/plotly makes it impossible to select the controls directly (drug, moa, targets)
# As a result, I remove this option entirely
controls = c("DMSO", "Bortezomib", "Irinotecan / SN-38")
selectable_drugs = tsne_dat %>% filter(isActive == "Active") %>% filter(!Drug %in% controls) %>%
  dplyr::pull(Drug) %>% unique() %>% sort()
selectable_moas = tsne_dat %>% filter(isActive == "Active") %>% filter(!Drug %in% controls) %>%
  dplyr::pull(MOA) %>% unique() %>% sort()
selectable_targets = tsne_dat %>% filter(isActive == "Active") %>% filter(!Drug %in% controls) %>%
  dplyr::pull(Target) %>% strsplit(c(",|;")) %>% unlist() %>% trimws() %>% toupper() %>% unique() %>% 
  sort()
selectable_lines = tsne_dat %>% filter(isActive == "Active") %>% filter(!Drug %in% controls) %>%
  dplyr::pull(Line) %>% unique() %>% sort()

xmin = min(tsne_dat$X)
xmax = max(tsne_dat$X)
ymin = min(tsne_dat$Y)
ymax = max(tsne_dat$Y)

color_scale = setNames(
  object = c("#b15928", "#1f78b4", "#fb9a99", "#33a02c"), 
  nm = sort(unique(tsne_dat$Line)))

shiny_theme <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(1, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")))}

server <- function(input, output, session) {
  output$tsneplot <- renderPlotly({
    drugs = input$selectedDrugs
    moas = input$selectedMOAs
    targets = input$selectedTargets
    lines = input$selectedLines
    
    unselected = tsne_dat
    selected = tsne_dat[tsne_dat$isActive == "Active", ]
    if(length(drugs) > 0) {
      drugs = unique(c(drugs, "DMSO", "Bortezomib", "Irinotecan / SN-38"))
      selected = selected[selected$Drug %in% drugs, ]
    }
    if(length(moas) > 0) selected = selected[selected$MOA %in% moas, ]
    if(length(targets) > 0) {
      sel_indices = targets %>% sapply(grep, x = toupper(selected$Target), fixed = TRUE) %>% 
        unlist() %>% unname() %>% unique()
      selected = selected[sel_indices, ]
    }
    if(length(lines) > 0) selected = selected[selected$Line %in% lines, ]

    dmso_control = unselected[unselected$Drug == "DMSO", ]
    
    gp = ggplot(mapping = aes(x = X, y = Y, text = sprintf(
      "Line: %s<br>Drug: %s<br>Target: %s<br>Mode of Action: %s", Line, Drug, Target, MOA))) + 
      geom_point(data = unselected, alpha = 0.5, color = "gray") + 
      geom_point(data = selected, mapping = aes(color = Line)) + 
      geom_point(data = dmso_control, mapping = aes(color = Line), size = 3, shape = 15) + 
      shiny_theme() + scale_color_manual(values = color_scale) + 
      theme(legend.position = "none") + labs(color = "", shape = "") + 
      coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
    gp = ggplotly(gp, tooltip = "text")
    gp$elementId = NULL
    gp
  })
}

ui <- fluidPage(
  titlePanel("Drug Effect Phenotypes"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("selectedDrugs", "Drugs",
                  c("Choose Drugs:" = "", selectable_drugs),
                  multiple = TRUE),
      selectInput("selectedMOAs", "MOAs",
                  c("Choose Mode of Action:" = "", selectable_moas),
                  multiple = TRUE),
      selectInput("selectedTargets", "Targets",
                  c("Choose Targets:" = "", selectable_targets),
                  multiple = TRUE),
      selectInput("selectedLines", "Lines",
                  c("Choose Lines:" = "", selectable_lines),
                  multiple = TRUE)
    ),
    mainPanel(plotlyOutput("tsneplot", height = WINDOW_SIZE, width = WINDOW_SIZE))
  )
)

shinyApp(ui = ui, server = server)
