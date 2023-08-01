rm(list=ls(all=TRUE))

library(shiny)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(UpSetR)



#Step 1: Load .txt files
#GWAS Information Table
setwd("~/Documents/PDGWAS")
Info_Table <- read.table('GWASInfoTable.txt',
                         header = TRUE,
                         sep = "\t",
                         dec = ".")


#all cell types
setwd("~/Documents/PDGWAS")
all_cell <- read.table('All_CellTypes.txt',              # TXT data file indicated as string or full path to the file
                       header = TRUE,    # Whether to display the header (TRUE) or not (FALSE)
                       sep = "",          # Separator of the columns of the file
                       dec = ".")  

#all genomic risk loci
setwd("~/Documents/PDGWAS")
all_genomic_loci <- read.table('All_Genomic_Loci.txt',              # TXT data file indicated as string or full path to the file
                       header = TRUE,    # Whether to display the header (TRUE) or not (FALSE)
                       sep = "",          # Separator of the columns of the file
                       dec = ".")  


#all independent sig snps
setwd("~/Documents/PDGWAS")
all_ind_snps <- read.table('AllIndSNPs.txt',              # TXT data file indicated as string or full path to the file
                           header = TRUE,    # Whether to display the header (TRUE) or not (FALSE)
                           sep = "",          # Separator of the columns of the file
                           dec = ".") 


#all snps
setwd("~/Documents/PDGWAS")
all_snps <- read.table('All_SNPs.txt',              # TXT data file indicated as string or full path to the file
                           header = TRUE,    # Whether to display the header (TRUE) or not (FALSE)
                           sep = "",          # Separator of the columns of the file
                           dec = ".") 







#Step 2: Reshape dataframes:

#####Info Table GWAS
#Reshape Info Table from a wide dataframe to a long dataframe
rownames_info <- as.list(names(Info_Table))
Info_Table <- t(Info_Table)
Info_Table <- as.data.frame(Info_Table)
rownames(Info_Table) <- rownames_info
names(Info_Table) <- Info_Table[1, ]
Info_Table <- Info_Table[-1, ]
Info_Table$GWAS <- rownames(Info_Table)
Info_Table<-Info_Table[,c(7, 1:6)]






#####Genomic Risk Loci
#Summarize by GWAS
upset_genomicloci <- all_genomic_loci %>%
  group_by(rsID, GWAS) %>%
  summarize(count = n()) %>%
  ungroup()

#Dummy Coding the Summarized Cell Types into Separate GWAS columns 
upset_genomicloci <- upset_genomicloci %>%
  pivot_wider(names_from = GWAS, values_from = count, values_fill = 0)

#Convert the tibble from pivot_wider function into a dataframe
wide_df <- as.data.frame(upset_genomicloci)

#Rename it back to Genomic Loci
upset_genomicloci <- wide_df




######Nearest Genes
#Summarize by GWAS
all_snps <- all_snps %>%
  arrange(desc(-gwasP))

sig_snps <- all_snps[all_snps$gwasP < (5*10^-8),]
sig_snps <- sig_snps[rowSums(is.na(sig_snps)) == 0, ]

upset_nearGene <- sig_snps %>%
  group_by(nearestGene, GWAS) %>%
  summarize(count = n()) %>%
  ungroup()

#Dummy Coding the Summarized Cell Types into Separate GWAS columns 
upset_nearGene <- upset_nearGene %>%
  pivot_wider(names_from = GWAS, values_from = count, values_fill = 0)

#Convert the tibble from pivot_wider function into a dataframe
wide_df <- as.data.frame(upset_nearGene)

#Extract Gene name column
nearGene_col <- wide_df$nearestGene

#Multiple of the same Cell Type within the same GWAS is recoded as 1 Cell Type
wide_df[wide_df > 1] = 1
wide_df$nearestGene <- nearGene_col

#Rename it back to Nearest Genes
upset_nearGene <- wide_df



###Significant Cell Types
#Subset all Cell Types with P.adj.pds < 0.05
all_cell$logP.adj.pds = -log10(all_cell$P.adj.pds)
all_cell <- all_cell %>%
  arrange(desc(logP.adj.pds))

celltypes_05 <- subset(all_cell, P.adj.pds <= 0.05)

#Summarize Cell Types by GWAS
celltypes_05 <- celltypes_05 %>%
  group_by(Cell_type, GWAS) %>%
  summarize(count = n()) %>%
  ungroup()

#Dummy Coding the Summarized Cell Types into Separate GWAS columns 
celltypes_05 <- celltypes_05 %>%
  pivot_wider(names_from = GWAS, values_from = count, values_fill = 0)

#Convert the tibble from pivot_wider function into a dataframe
wide_df <- as.data.frame(celltypes_05)

#Multiple of the same Cell Type within the same GWAS is recoded as 1 Cell Type
wide_df[wide_df > 1] = 1

#Rename it back to Cell Types of P<0.05
celltypes_05 <- wide_df




#Step 3: Reshape cell types dataset for the barplot
data <- all_cell
data_colsubset = data %>%
  select("Dataset", "Cell_type", "P.adj.pds", "GWAS")

data_colsubset$logP.adj.pds = -log10(data_colsubset$P.adj.pds)
data_colsubset <- data_colsubset %>%
  arrange(desc(logP.adj.pds))

n_row_for_plot = 50
data_colsubset_rowsubset = data_colsubset[1:n_row_for_plot, ]

data_colsubset_rowsubset$Cell_type_GWAS = paste(data_colsubset_rowsubset$Cell_type, data_colsubset_rowsubset$GWAS, data_colsubset_rowsubset$Dataset, sep="_")

data_colsubset_rowsubset$Cell_type_GWAS <- factor(data_colsubset_rowsubset$Cell_type_GWAS, levels = data_colsubset_rowsubset$Cell_type_GWAS[order(data_colsubset_rowsubset$logP.adj.pds)])

head(data_colsubset)
head(data_colsubset_rowsubset)




#Step 4: Create lists containing GWAS names for Upset plots
gwas_names <- c(names(celltypes_05)[2:5])
print(gwas_names)
#[1] "Simon"        "Nalls"        "Pankratz"     "Blauwendraat"




#Step 5: Run our Shiny app to visualize GWAS comparison '
######################
#RShiny Tool: Cell Type Table, Histogram of P.adj.pds-values, Barplot of P.adj.pds-values, and Upset Plot of SNP Overlap
######################
# User Interface: 

ui <- fluidPage(
  div(  # Header Panel with blue color and text
    class = "header-panel",
    style = "background-color: #302a75; color: #FFFFFF; padding: 10px;",
    "GWAS annotation tools can be found on the FUMA platform at fuma.ctglab.nl. Thanks FUMA!"
  ),  # Add the header panel
  titlePanel("Explorer Tool for Parkinson's Disease GWAS"), #Title Panel
  sidebarLayout(
    sidebarPanel(  #5 inputs: 1) Cell Types Table Column Input, 2) GWAS Input, 3) scRNAseq Dataset Input, 4) Number of rows to display in Cell Type Table, 5) Barplot Height
      selectInput("col_input", "1. Select Columns to Display in Cell Types Table:", choices = c("All", names(all_cell)), multiple = TRUE),
      selectInput("gwas_input", "2. Select GWAS:", choices = c("All", unique(all_cell$GWAS))),
      selectInput("dataset_input", "3. Select scRNAseq Dataset:", choices = c("All", unique(all_cell$Dataset))),
      numericInput("row_input", "4. Select Number of Cell Types to Display:", value = 10, min = 1),
      numericInput("barplot_height", "5. Adjust Barplot Height:", value = 400, min = 100)
    ),
    mainPanel( #7 Output Panels: 1) Cell Types Table, 2) Cell Types Barplot, 3) Simple Manhattan Plot 4) GWAS Information Table, 5) Lead SNPs Upset, 6) Ind Sig SNPs Upset, 7) Cell Types Upset plot
      tabsetPanel(
        tabPanel("Cell Types Table", 
                 div(
                   style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   div(
                     style = "margin-bottom: 20px; text-align: center;",
                     "Cell Types Table shows the cell types with lowest p-values (adjusted within dataset). Select the columns you would like to see in Input Box #1 along with the column Cell Type, which shows you cell types, and the column P.adj.pds, which shows you the p-values adjusted within the cell type's scRNAseq dataset. You can also select the column GWAS to see which GWAS signal is associated with this cell type, and the column Dataset to see which scRNAseq dataset contains this cell type. Other columns include BETA (beta coefficient), BETA_STD (standardized beta coefficient), SE (standard error), P (unadjusted p-value), and P.adj (P-value adjusted between all scRNAseq datasets). You can also select GWAS in Input Box #2 and scRNAseq dataset in Input Box #3 to look at only the cell types from these groups. The number of cell types displayed in this table can be adjusted in Input Box #4." 
                   ),
                   dataTableOutput("table"),
                   style = "margin-top: 20px;"  # Add some spacing between the text and the barplot
                 )
        ),
        tabPanel("Cell Types Barplot",
                 div(
                   style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   div(
                     style = "margin-bottom: 20px; text-align: center;",
                     "Cell Types Barplot shows you how many cell types are considered significant and which cell types those are. Cell Types are listed on the y-axis of this plot and -log10 transformed P.adj.pds values (adjusted within scRNAseq dataset) on the x-axis. Cell types that fall past the significance threshold (blue dotted line) of -log10(5x10^-8) are considered significant. You can select cell types identified from a specific GWAS in Input Box #2 and cell types from a specific scRNAseq dataset in Input Box #3. Make sure to adjust the height of this barplot in Input Box #5 to be able to see the y-axis as the plot changes."
                   ),
                   plotOutput("barplot", height = "400px")
                 )
        ),
        tabPanel("Nearest Genes Upset Plot",
                 div(
                   style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   div(
                     style = "margin-bottom: 20px; text-align: center;",
                     "Nearest Genes Upset Plot displays the overlap in all the genes nearest to significant SNPs identified from each Parkinson's Disease GWAS, with GWAS displayed on the x-axis and the number of nearest genes displayed on the y-axis. Nearest Genes that overlap are listed below the Upset plot. Nearest genes were identified using SNP2GENE on the FUMA platform."
                   ),
                   plotOutput("upsetNearGenePlot"),
                   tableOutput("NearGeneTable")
                 )
        ),
        tabPanel("Genomic Risk Loci Upset Plot",
                 div(
                   style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   div(
                     style = "margin-bottom: 20px; text-align: center;",
                     "Genomic Risk Loci Upset plot displays the overlap in genomic risk loci identified from each Parkinson's Disease GWAS, with GWAS displayed on the x-axis and number of genomic risk loci displayed on the y-axis. Genomic rick loci that overlap are listed below the Upset plot. Genomic risk loci were identified using SNP2GENE on the FUMA platform; SNP2GENE lists genomic risk loci by the rsID of Lead SNP found within the risk loci. SNP2GENE defines Lead SNPs as SNPs that are significant at the genome-wide threshold level of 5x10^-8, and which are independent of each other by being below a correlation (r^2) threshold of 0.1."
                   ),
                   plotOutput("upsetGenomicLociPlot"),
                   tableOutput("GenomicLociTable")
                 )
        ),
        tabPanel("Cell Types Upset Plot",
                 div(
                   style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;",
                   div(
                     style = "margin-bottom: 20px; text-align: center;",
                     "Cell Types Upset Plot displays whether there is an overlap in significant cell types identified from each Parkinson's Disease GWAS, with GWAS displayed on the x-axis and cell types displayed on the y-axis. Cell types that overlap are listed below the Upset plot. Cell types were identified using Cell Type on the FUMA platform; Cell Type analyzes genes that emerge from GWAS signal and uses MAGMA gene-property analysis to correlate gene expression profiles with cell types. Each cell type's gene expression profile is provided from single-cell RNA sequencing datasets. All cell types displayed fall below a P.adj.pds (adjusted within scRNAseq dataset) threshold of 0.05."
                   ),
                   plotOutput("upsetPlot"),
                   tableOutput("celltypesTable")
                 )
        ),
        tabPanel("GWAS Information Table", tableOutput("infoTable")),  # New tab panel for Info Table
      )
    )
  )
)

# Server part

server <- function(input, output) {
  
  filtered_data <- reactive({ #filtered reaction for all_cell dataframe which contains cell types, allows subsetting by GWAS and dataset
    data <- all_cell
    if (input$gwas_input != "All") {
      data <- data[data$GWAS == input$gwas_input,]
    }
    if (input$dataset_input != "All") {
      data <- data[data$Dataset == input$dataset_input,]
    }
    data <- data[order(data$P.adj.pds),]
    data[1:input$row_input, input$col_input]
  })
  
  output$table <- renderDataTable({ #Outputs the 1) Cell Types Table
    filtered_data()
  })
  
  output$scatterplot <- renderPlot({ #Outputs the 3) Simple Manhattan Plot
    data <- all_snps
    
    if (input$gwas_input != "All") {
      data <- data[data$GWAS == input$gwas_input, ]
    }
    
    plot(data$pos, data$gwasP,
         xlab = "Position (pos)",
         ylab = "GWAS P-value (gwasP)",
         main = "Scatterplot of Position (pos) vs GWAS P-value (gwasP)")
    
    abline(h = 5*(10^-8), col = 'blue', lwd = 3, lty = 'dashed')
  })
  
  output$barplot <- renderPlot({  #Outputs the 2) Cell Types Barplot
    data <- data_colsubset_rowsubset
    x_col <- "logP.adj.pds"
    y_col <- "Cell_type_GWAS"
    
    if (input$gwas_input != "All") {
      data <- data[data$GWAS == input$gwas_input, ]
    }
    
    data <- data[order(data$logP.adj.pds), ]  
    data <- unique(data)
    data[, y_col] <- factor(data[, y_col], levels = data[, y_col][order(data$logP.adj.pds)])  
    
    ggplot(data, aes(x = get(x_col), y = get(y_col))) +
      geom_bar(stat = "identity") +
      labs(title = "Barplot of Cell Type Associations with PD",
           x = "-log10(P.adj.pds)",
           y = "Cell Type") +
      geom_vline(xintercept = -log10(0.05), col = "blue")
  }, height = function() {
    validate(need(!is.null(input$barplot_height), "Please provide a barplot height."))
    input$barplot_height
  })
  
  
  output$upsetPlot <- renderPlot({    # Outputs the 7) Cell Types Upset
    upset(celltypes_05, sets = gwas_names,
          keep.order=TRUE,
          sets.bar.color = "#56B4E9", order.by = "freq",
          matrix.color = "black", main.bar.color = "black",
          mainbar.y.label = "Number of Significant Cell Types")
  })
  
  output$celltypesTable <- renderTable({ # Outputs the overlapping cell types in 7)
    # Subset the data frame to include only numeric columns
    cell_types_numeric <- celltypes_05[, sapply(celltypes_05, is.numeric)]
    
    # Calculate the row sums for the numeric columns
    celltypes_05$sum <- rowSums(cell_types_numeric)
    celltypes_05many <- celltypes_05[celltypes_05$sum > 1, ]
    head(celltypes_05many)
  })
  
  output$upsetNearGenePlot <- renderPlot({    # Outputs the 4) nearest Gene UpSet plot
    upset(upset_nearGene, sets = gwas_names,
          keep.order=TRUE,
          sets.bar.color = "#56B4E9", order.by = "freq",
          matrix.color = "black", main.bar.color = "black",
          mainbar.y.label = "Number of Genes near Genome-Wide Significant SNPs")
  })
  
  output$NearGeneTable <- renderTable({ # Outputs the overlapping near genes in 4)
    # Subset the data frame to include only numeric columns
    near_gene_numeric <- upset_nearGene[, sapply(upset_nearGene, is.numeric)]
    # Calculate the row sums for the numeric columns
    upset_nearGene$sum <- rowSums(near_gene_numeric)
    
    upset_nearGene_many <- upset_nearGene[upset_nearGene$sum > 1, ]
    upset_nearGene_many <- upset_nearGene_many %>%
      arrange(desc(sum))
    print(upset_nearGene_many)
  })
  
  output$upsetGenomicLociPlot <- renderPlot({    # Outputs the 5) Genomic Loci UpSet plot
    upset(upset_genomicloci, sets = gwas_names,
          keep.order=TRUE,
          sets.bar.color = "#56B4E9", order.by = "freq",
          matrix.color = "black", main.bar.color = "black",
          mainbar.y.label = "Number of Genomic Risk Loci")
  })
  
  output$GenomicLociTable <- renderTable({ # Outputs the overlapping Genomic Loci in 5)
    # Subset the data frame to include only numeric columns
    gen_loci_numeric <- upset_genomicloci[, sapply(upset_genomicloci, is.numeric)]
    # Calculate the row sums for the numeric columns
    upset_genomicloci$sum <- rowSums(gen_loci_numeric)
    
    upset_genomicloci_many <- upset_genomicloci[upset_genomicloci$sum > 1, ]
    print(upset_genomicloci_many)
  })
  
  output$infoTable <- renderTable({ #Outputs "3) GWAS Information Table"
    Info_Table
  })
}

shinyApp(ui = ui, server = server)









