.cleanFile <- function (file)
{
    file.vector <- unlist(strsplit(file, "\\."))
    if (file.vector[length(file.vector)] %in% c("txt",
        "sif", "tab", "XGMML",
        "NOA", "EDA", "xgmml", "eda",
        "noa", "net", "pdf")) {
        file.vector <- file.vector[-length(file.vector)]
        file <- paste(file.vector, collapse = ".")
    }
    return(file)
}
