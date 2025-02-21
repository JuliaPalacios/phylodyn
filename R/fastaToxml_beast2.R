#' Generate xml files from fasta files
#' 
#' @description Given a fasta file from seqgen(), this function outputs an xml file as if we run BEAUti
#' @param inputfasta The path to the input fasta file
#' @param outputxml The path to the output xml file
#' @param xml_template The path to the xml template
#' @param clock_rate The clock rate used in BEAUti
#' @return An xml file
#' @examples 
#' fastaToxml_beast2("/path/inputfasta.fasta", "/path/outputxml.xml", "/path/xml_template.txt", clock_rate = 2) 
#' @export
fastaToxml_beast2 <- function(inputfasta, outputxml, xml_template, clock_rate) {
  library(seqinr)
  
  clock_rate <- sprintf("%0.1f", clock_rate)
  line_clock_rate <- 62 # this is fixed
  
  # check if the output file exist. Delete file if it exists.
  if (file.exists(outputxml)) {
    file.remove(outputxml)
  }
  
  # read template.txt
  template_open <- file(xml_template, open="r")
  template_lines <- readLines(template_open) 
  lines_total <- length(template_lines)
  lines_before <- 6  # This is fixed. the number of lines before the sequence data
  
  # write template before data
  for (i in 1:lines_before){
    cat(paste0(template_lines[i], "\n"), file = outputxml, append = T)
  }
  
  # fasta file data
  alig <- read.fasta(inputfasta, as.string = T, forceDNAtolower = F)
  n <- length(alig)
  ids <- names(alig)
  
  for(i in seq(1, n)){
    id <- ids[i]
    seq_id <- paste0("seq_", id)
    sequence <- alig[[i]][1]
    line <- sprintf('\t\t<sequence id="%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n', seq_id, id, sequence)
    cat(line, file = outputxml, append = T)
  }
  
  # write template after data
  for (i in (lines_before+1):lines_total){
    if(i == line_clock_rate){  # for different clock rate
      line <- sprintf(paste0(template_lines[62], "\n"), clock_rate)
      cat(line, file = outputxml, append = T)
    }else{
      cat(paste0(template_lines[i], "\n"), file = outputxml, append = T)
    }
  }
  
  close(template_open)
}








