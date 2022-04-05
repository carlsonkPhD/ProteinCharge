#' Calculate theoretical charge
#'
#' Calculate the theoretical charge of the given protein sequence over a range of pH values
#' @param sequence The single letter abbreviated protein sequence. A string or character vector of length one.
#' @param read_source File path to the text file containing the protein sequence.
#' @param state whether to calculate the charge under oxidizing or reducing conditions. default is reducing.
#' @param chart_title The title to be displayed by ggplot2 above the produced graph
#' @param data_type Instructs function whether to consider the sequence or read_source arguments. default value is "string".
#' @param pH_start Starting point for theoretical titration. Must be numeric.
#' @param pH_stop End point of theoretical titration. Must be numeric
#' @return Returns the protein charge to Stdout and creates graph in Plots panel
#' titration(sequence = "DYKDDDDK", chart_title = "Flag peptide")
#' titration(read_source = "flagPeptide.txt", chart_title= "Flag peptide", data_type = "file")
#' @import ggplot2
#' @import stringr
#' @export
titration <- function(sequence = "", read_source = "", state = c("reduced", "oxidized"), chart_title = "", data_type = c("string", "file"), pH_start=4.0, pH_stop=10.0){
     data_type <- match.arg(data_type)
     state <- match.arg(state)
     
     if (data_type == "string") {
          vector_aasequence <- sequence
     } else if (data_type == "file") {
          vector_aasequence <- as.character(unlist(read.table(read_source)))
     }
     
     if (pH_stop - pH_start <= 0) {
          cat("Please ensure the pH_start value is creater than pH_stop value\n")
          opt <- options(show.error.messages = FALSE)
          on.exit(options(opt))
          stop()
     }
     
     
     num_Aspartic <- str_count(vector_aasequence, "D")
     num_Glutamic <- str_count(vector_aasequence, "E")
     num_Histidine <- str_count(vector_aasequence, "H")
     num_Cysteine <- str_count(vector_aasequence, "C")
     num_Tyrosine <- str_count(vector_aasequence, "Y")
     num_Lysine <- str_count(vector_aasequence, "K")
     num_Arginine <- str_count(vector_aasequence, "R")
     
     pH_seq <- seq(from = pH_start, to = pH_stop, by = 0.1)
     charges_vector <- c()
     
     for (i in pH_seq){
          #EMBOSS values
          aspartic_charge <- 1/(1 + 10^(3.9-i))
          glutamic_charge <- 1/(1 +10^(4.1-i))
          cysteinie_charge <- 1/(1 +10^(8.5-i))
          tyrosine_charge <- 1/(1 +10^(10.1-i))
          cterm_charge <- 1/(1 +10^(3.6-i))
          
          histidine_charge <- 1/(1 +10^(i-6.5))
          lysine_charge <- 1/(1 +10^(i-10.8))
          arginine_charge <- 1/(1 +10^(i-12.5))
          nterm_charge <- 1/(1 +10^(i-8.6))
          if (state == "reduced"){
               total_charge <- ((num_Histidine * histidine_charge) + (num_Lysine * lysine_charge) + 
                                     (num_Arginine * arginine_charge) + nterm_charge) - ((num_Aspartic * aspartic_charge) 
                                                                                         + (num_Glutamic * glutamic_charge) + (num_Cysteine * cysteinie_charge) 
                                                                                         + (num_Tyrosine * tyrosine_charge) + cterm_charge)
          } else if (state == "oxidized") {
               total_charge <- ((num_Histidine * histidine_charge) + (num_Lysine * lysine_charge) + 
                                     (num_Arginine * arginine_charge) + nterm_charge) - ((num_Aspartic * aspartic_charge) 
                                                                                         + (num_Glutamic * glutamic_charge) + (num_Tyrosine * tyrosine_charge) + cterm_charge)
          }
          
          charges_vector <- c(charges_vector, total_charge)
     }
     my_df <- data.frame(pH = pH_seq, charge = charges_vector)
     print(my_df)
     
     p <- ggplot(my_df, aes(x = pH, y = charge)) + geom_point() + labs(title = chart_title)
     p.theme <- p + theme(panel.grid.major.y = element_line(color = "blue"))
     p.scale <- p.theme + scale_x_continuous(breaks = c(pH_start:pH_stop))
     p.labels <- p.scale + labs(x = "pH", y = "Charge")
     p.labels
}