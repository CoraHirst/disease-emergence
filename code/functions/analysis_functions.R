# function to pull each sheet of an excel file as its own data frame and store all DFs in a list
multiplesheets <- function(excel.file) { 
  
  #load necessary packages
  require(readxl)
  library(readxl)
  
  # getting info about all excel sheets 
  sheets <- excel_sheets(excel.file) 
  tibble <- lapply(sheets, function(x) read_excel(excel.file, sheet = x)) 
  data_frames <- lapply(tibble, as.data.frame) 
  
  # assigning names to data frames 
  names(data_frames) <- sheets 
  
  return(data_frames)
} 