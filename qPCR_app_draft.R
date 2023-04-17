library(gtsummary)
library(tidyverse)
library(janitor)
library(gridExtra)
library(readxl)
library(ggpmisc)
library(tidyquant)
library(patchwork)
library(forestmangr)
library(kableExtra)
library(gt)
library(scales)
library(dplyr)
library(ggplot2)
library(shiny)
library(shinythemes)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("flatly"),
  navbarPage(
    "qPCR",
    tabPanel("qPCR 1",
             sidebarPanel(
               fileInput("file", "Choose Excel file"),
               selectInput("sheet", "Select sheet", choices = NULL),
               actionButton("process", "Process Data"),
               actionButton("overview", "Overview of Data"),
               selectInput("tissue", "Select tissue", choices = NULL),
               actionButton("b_t", "Tissue sample"),
               actionButton("b_p", "Plasma sample")
               
               
             ),
             mainPanel(
               tags$style(type="text/css",
                          ".shiny-output-error { visibility: hidden; }",
                          ".shiny-output-error:before { visibility: hidden; }"),
               plotOutput("cq_conc_plot"),
               tableOutput("biod_tissues_table"),
               plotOutput("biod_tissues_plot"),
               tableOutput("biod_plasma_table"),
               plotOutput("biod_plasma_plot")
             )
             
    ),
    tabPanel("Bar 2", "This is intentionally left blank"),
    tabPanel("Bar 3", "This is intentionally left blank")
   )
 )



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(scipen = 100, digits = 4)
  theme_set(theme_tq(base_size = 14))
  
  # Read Excel data
  read_excel_data <- function(file_name, sheet_name) {
    df <- readxl::read_xlsx(file_name, sheet = sheet_name)
    df <- df %>% row_to_names(row_number = 23) %>% clean_names()
    return(df)
  }
  
  # Preprocess data
  preprocess_data <- function(df){
    processed_data <- df %>%
      mutate(cq = as.numeric(cq), quantity = as.numeric(quantity)) %>%
      drop_na(cq, quantity) %>%
      select(3, 4, 11) %>% 
      separate(col = sample, into = c("tissue", "type", "quantity", "rep"), sep = "-") %>% 
      mutate(tissue1 = ifelse(tissue == "Bra", "Brain",
                              ifelse(tissue == "Hea", "Heart",
                                     ifelse(tissue == "Kid", "Kidney",
                                            ifelse(tissue == "Lun", "Lungs", 
                                                   ifelse(tissue == "Pla", "Plasma",
                                                          ifelse(tissue == "SkM", "Muscle",
                                                                 ifelse(tissue == "Spl", "Spleen", "Blood Cells")))))))) %>%
      mutate(quantity = as.numeric(quantity)) %>%
      ungroup()
    return(processed_data) 
  }
  
  #Cq vs conc
  cq_conc <- function(processed_data){
    s1 <- processed_data %>% filter(type=="Std") %>% group_by(quantity,tissue1)  %>% summarise(m=mean(cq)) %>%
      rename(input=quantity)  %>% mutate(lconc=log(input)) 
    
    #Single tissue line graph
    g1 <- ggplot(s1,aes(log(input),m)) + geom_point(size=2) + geom_smooth(method = "lm",se=F) + labs(x="Log Concentration     (ng/mL)", y="Cq")+theme_bw() + stat_poly_eq(coef.digits = 4,label.x = c(0.9,0.9,0.9),label.y = c(0.9,0.7,0.6)) + facet_wrap(~tissue1)
    
    #Line graph for all tissues
    g2 <- ggplot(s1,aes(log(input),m,color=tissue1,shape=tissue1)) + geom_point(size=2) + geom_smooth(method = "lm",se=F) + labs(x="Log Concentration (ng/mL)", y="Cq")+theme_bw()
    
    return(g1 / g2)
  }
  
  
  #Bio distribution in Tissues
  biod_tissues <- function(tissue){
    
    # Get processed data
    data <- processed_data()
    
    # Filter data based on selected tissue
    s1 <- data %>% filter(tissue=={{tissue}} & type=="Std") %>% group_by(quantity,tissue1)  %>% summarise(m=mean(cq)) %>% rename(input=quantity)  %>% mutate(lconc=log(input))
    a <- unique(s1$tissue1)
    
    #Max cq
    maxcq <- 27
    
    mod <- lm(lconc~m,data = s1)
    
    s1$bconc0 <- predict(mod,newdata = s1)
    s1 <- s1 %>% mutate(bconc=exp(bconc0),pdiff=(bconc-input)/input*100) %>% select(-bconc0) %>% round_df(digits = 3,rf = "signif")
    
    colnames(s1) <- c("Input (ng/mL)","Tissue","Mean Cq","Log Input (ng/mL)","Back Calculated \n Concentration (ng/mL)","Percent \n Difference")
    
    desc_table <- s1 %>% gt()
    
    sam1 <- data %>% filter(tissue=={{tissue}} & type=="Ukn") %>% rename(sample=quantity,m=cq) %>% filter(m<maxcq)
    sam1$conc <- exp(predict(mod,newdata = sam1))
    sam_abs <- sam1  %>% mutate(conc_nm=conc*100,conc_ug_g=conc_nm*445.5/1000) #mol.wt = 445500Da; tissue conc = 10mg/mL
    
    sam.s <- sam_abs %>% group_by(sample) %>% summarise(m=mean(conc_ug_g),sd=sd(conc_ug_g)) %>% mutate(sample=as.character(sample))
    
    sam.ss <- sam.s %>% ungroup() %>%  summarise(im=mean(m),sd=sd(m),cv=sd/im*100) %>% mutate(Tissue=paste0(a))
    
    overall <- data.frame(sample="Mean (±SD)",m=sam.ss$im,sd=sam.ss$sd)
    overall1 <- overall %>% mutate(sample=paste0(a))
    
    sam.s$sd <- 0
    sam.sss <- bind_rows(sam.s,overall)
    
    # Plotting the graph
    ggplot(sam.sss,aes(sample,m)) + geom_errorbar(aes(ymin=m,ymax=m+sd),width=0.2,color="blue",alpha=0.6) + geom_col(col="blue",alpha=0.7,fill="blue") + labs(y="Concentration (µg/g)",x="Animal ID",title = paste0("Tissue=",a),caption = "Dose: 1 mg/kg") +
      scale_y_continuous(breaks = pretty_breaks(n = 5))
    
    
    
    
  }
  
  biod_plasma <- function(tissue){
    
    s1 <- processed_data %>% filter(tissue=={{tissue}} & type=="Std") %>% group_by(quantity,tissue1)  %>%
      summarise(m=mean(cq)) %>% rename(input=quantity)  %>% mutate(lconc=log(input)) 
    a <- unique(s1$tissue1)
    maxcq <- 27
    
    mod <- lm(lconc~m,data = s1)
    
    #Predicting the bconc 0
    s1$bconc0 <- predict(mod,newdata = s1)
    
    #
    s1 <- s1 %>% mutate(bconc=exp(bconc0),pdiff=(bconc-input)/input*100) %>% select(-bconc0) %>% round_df(digits = 3,rf = "signif")
    
    colnames(s1) <- c("Input (ng/mL)","Tissue","Mean Cq","Log Input (ng/mL)","Back Calculated \n Concentration (ng/mL)","Percent \n Difference")
    
    
    desc_table <- s1 %>% gt()
    sam1 <- processed_data %>% filter(tissue=={{tissue}} & type=="Ukn") %>% rename(sample=quantity,m=cq) %>% filter(m<maxcq)
    sam1$conc <- exp(predict(mod,newdata = sam1))
    sam_abs <- sam1  %>% mutate(conc_nm=conc,conc_ug_g=conc_nm*445.5/1000)
    
    sam.s <- sam_abs %>% group_by(sample) %>% summarise(m=mean(conc_ug_g),sd=sd(conc_ug_g)) %>%
      mutate(sample=as.character(sample))
    
    
    
    
    sam.ss <- sam.s %>% ungroup() %>%  summarise(im=mean(m),sd=sd(m),cv=sd/im*100) %>% mutate(Tissue=paste0(a))
    
    overall <- data.frame(sample="Mean (±SD)",m=sam.ss$im,sd=sam.ss$sd)
    
    overall1 <- overall %>% mutate(sample=paste0(a))
    
    sam.s$sd <- 0
    sam.sss <- bind_rows(sam.s,overall) %>% mutate(m=m*1000,sd=sd*1000)
    
    
    ggplot(sam.sss,aes(sample,m)) + geom_errorbar(aes(ymin=m,ymax=m+sd),width=0.2,color="blue",alpha=0.6) +
      geom_col(col="blue",alpha=0.7,fill="blue") + labs(y="Concentration (ng/mL)",x="Animal ID",title = paste0("Matrix=",a),caption = "Dose: 1 mg/kg") + scale_y_continuous(breaks = pretty_breaks(n = 5)) 
    
    
  }
  
  #Update the select input choices for sheet based on uploaded file
  observe({
    req(input$file)
    sheets <- readxl::excel_sheets(input$file$datapath)
    if (is.null(sheets)) {
      showNotification("Error: No sheets found in uploaded file", type = "error")
    } else {
      updateSelectInput(session, "sheet", choices = sheets)
    }
  })
  
  
  #Process the data when the process button is clicked
  processed_data <- eventReactive(input$process, {
    req(input$file, input$sheet)
    df <- readxl::read_xlsx(input$file$datapath, sheet = input$sheet)
    df <- df %>% row_to_names(row_number = 23) %>% clean_names()
    processed_data <- df %>%
      mutate(cq = as.numeric(cq), quantity = as.numeric(quantity)) %>%
      drop_na(cq, quantity) %>%
      select(3, 4, 11) %>%
      separate(col = sample, into = c("tissue", "type", "quantity", "rep"), sep = "-") %>%
      mutate(tissue1 = case_when(
        tissue == "Bra" ~ "Bra",
        tissue == "Hea" ~ "Hea",
        tissue == "Kid" ~ "Kid",
        tissue == "Lun" ~ "Lun",
        tissue == "Pla" ~ "Pla",
        tissue == "SkM" ~ "SkM",
        tissue == "Spl" ~ "Spl",
        TRUE ~ "Blood Cells"
      )) %>%
      mutate(quantity = as.numeric(quantity)) %>%
      ungroup()
    processed_data
  })
  
  #Update the select input choices for tissue based on processed data
  observe({
    req(processed_data())
    tissues <- processed_data() %>%
      select(tissue1) %>%
      pull()
    updateSelectInput(session, "tissue", choices = unique(tissues))
  })
  
  #Render the Cq vs conc plot
  # Use observeEvent to trigger plot generation when overview button is clicked
  observeEvent(input$overview, {
    output$cq_conc_plot <- renderPlot({
      processed_data() %>% cq_conc()  
    })
  })
  
  # Update tissue plot when "Tissue sample" button is clicked
  observeEvent(input$b_t, {
    output$biod_tissues_plot <- renderPlot({
      biod_tissues(input$tissue)
    })
    
    output$biod_tissues_table <- render_gt({
      biod_tissues(input$tissue)$desc_table
    })
  })
  
  #For the Plasma
  observeEvent(input$b_p, {
    output$biod_tissues_plot <- renderPlot({
      biod_tissues(input$tissue)
    })
    
    output$biod_tissues_table <- render_gt({
      biod_tissues(input$tissue)$desc_table
    })
  })
    

    
}

# Run the application 
shinyApp(ui = ui, server = server)
