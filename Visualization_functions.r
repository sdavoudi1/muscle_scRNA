# This file contains all the functions that were used to create the images for the 
# scRNA Seq project.

# ----------------------------------------------------------------------------------------------

# This function creates a heatmap in which the values for each row are normalized between 0 and 1

# NEEDS WORK. DO NOT USE.
# 

createHeatmap <- function (input, max_color = "red", min_color = "peachpuff", show_scale = F){


    # Check if the package installed
    if (!require("plotly")) {install.packages("plotly"); require(plotly)}


    #Read in .csv info for heatmap
    nba <- read.csv(input, sep=",")


    #Make first column name of each row, take values from column [B Cells - Tendon] as table
    row.names(nba) <- nba$Line
    nba <- nba[,2:10]
    nba_matrix <- data.matrix(nba)


    #Sort matrix diagonally by the largest value of each row (top left to bottom right)
    row.max <- apply(nba_matrix,1,which.max)
    nba_matrix  <- nba_matrix[names(sort(row.max)),]


    # Load plotly and create graph
    library(plotly)
    var <- plot_ly(x=colnames(nba_matrix), y=rownames(nba_matrix), z = nba_matrix, type = "heatmap", colors = colorRamp(c(min_color, max_color)), showscale = show_scale)


    #Export plot using webshot then opens image
    if (!require("webshot")) install.packages("webshot")
    tmpFile <- tempfile(fileext = ".png")
    export(var, file = tmpFile)
    browseURL(tmpFile)

}
