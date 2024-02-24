#' m.shortcuts: A personal package to make my life easier
#'
#' The m.shortcuts package has a disparage collection of functions used in
#' common bioinformatics pipelines, and work from
#' [Miguel Angel Garcia-Campos](https://angelcampos.github.io/).
#'
#' @docType package
#' @name m.shortcuts
#' @aliases m.shortcuts-package
#' @keywords ggplot2
#' @importFrom rlang .data
"_PACKAGE"
NULL

#' Improved list of objects
#'
#' Prints a table showing the size of biggest objects in the global environment
#' (by default). Based on code by Dirk Eddelbuettel at
#' \href{https://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session}{StackOverflow}
#'
#' @param pos numeric. an alternative argument to name for specifying the environment as
#' a position in the search list. Mostly there for back compatibility.
#' @param pattern character. pattern to match objects listed
#' @param order.by character. Name of column to order results by. Set to NULL for not ordering
#' @param decreasing logical. Set to FALSE if increasing order is desired
#' @param head logical. By default limits results to n results
#' @param n numeric. Number of results to show when setting head to TRUE
#'
#' @return data.frame
#' @export
#'
#' @author Dirk Eddelbuettel
lsos <- function(pos = 1, pattern, order.by = "Size_Mb", decreasing = TRUE,
                         head = TRUE, n = 10){
    napply <- function(names, fn) sapply(names, function(x)
        fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- round(napply(names, utils::object.size) / 1024^2, digits = 2)
    obj.dim <- t(napply(names, function(x)
        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size_Mb", "Rows", "Columns")
    if(!is.null(order.by)){
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    }
    if (head)
        out <- head(out, n)
    out
}


#' Remove tmp objects from global environment
#'
#' @param collectGarbage logical. Activate garbage collector
#'
#' @return matrix By default outputs the gc() output matrix
#' @export
rmTmp <- function(collectGarbage = TRUE){
    tmpOL <- grep(x = ls(pos = ".GlobalEnv"), pattern = "^tmp", value = TRUE)
    rm(list = tmpOL, pos = ".GlobalEnv")
    if(collectGarbage){gc()}
}

#' Check figs
#'
#' Look for matches figures in both .R scripts and .Rmd notebooks
#'
#' @param figsDir Figures directory. Use empty
#' @param scripts character. File names of scripts to check
#' @param notebooks character. File names of Notebooks to check
#' @param fig_exts character. Extensions associated to figures, e.g. 'pdf',
#'  'jpg'. Do not include a "." before the extension.
#'
#' @return data.table
#' @export
check_figs <-  function(figsDir = "figs",
                        scripts = list.files(pattern = "\\.R$"),
                        notebooks = list.files(pattern = "\\.Rmd$"),
                        fig_exts = c("pdf", "png", "jpg")){
    if(!dir.exists(figsDir)){stop("Directory '", figsDir, "' does not exist")}
    figsInDir <- lapply(fig_exts, function(ext){
        tmp <- list.files(figsDir, ext, full.names = TRUE, recursive = TRUE)
        data.table::data.table(fig = tmp)
    }) %>% do.call(what = "rbind")
    figsScrpt <- lapply(scripts, function(file){
        script <- readLines(file)
        figFiles <- lapply(fig_exts, function(figExt){
            stringr::str_extract_all(script, pattern = file.path(figsDir, paste0("[[:graph:]]+\\.", figExt))) %>% unlist()
        }) %>% unlist() %>% unique()
        lapply(figFiles, function(figFile){
            data.table::data.table(fig = figFile, generated_by = file, generated_line = which(stringr::str_detect(script, figFile)))
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind")
    if(nrow(figsInDir) > 0){
        figsTable <- dplyr::full_join(figsScrpt, figsInDir, by = "fig")
    }else{
        figsTable <- figsScrpt
    }

    figs_dupGen <- figsTable[duplicated(dplyr::select(figsTable, c("fig", "generated_by"))),]$fig

    figsNoteB <- lapply(notebooks, function(file){
        notebook <- readLines(file)
        figFiles <- lapply(fig_exts, function(figExt){
            stringr::str_extract_all(notebook, pattern = file.path(figsDir, paste0("[[:graph:]]+\\.", figExt))) %>% unlist()
        }) %>% unlist() %>% unique()
        lapply(figFiles, function(figFile){
            if(any(stringr::str_detect(notebook, figFile))){
                data.frame(fig = figFile, used_in = file, used_line = stringr::str_which(notebook, figFile))
            }else{
                NULL
            }
        }) %>% do.call(what = "rbind")
    }) %>% do.call(what = "rbind")

    if(is.null(figsNoteB)){
        figsTable$used_in <- NA
        figsTable$used_line <- NA
    }else{
        figsTable <- dplyr::full_join(figsTable, figsNoteB, by = "fig")
    }
    figs_notUsed <- figsTable[is.na(figsTable$used_in),]$fig
    tmpT <- unique(dplyr::select(figsTable, c("fig", "used_in", "used_line")))
    figs_dupUsed <- tmpT[duplicated(dplyr::select(tmpT, c("fig", "used_in"))),]$fig
    if(length(figs_notUsed) > 0){
        warning("Found figures in dir, but not in code:\n",
                paste(figs_notUsed, collapse = "\n"))
    }
    if(length(figs_dupGen) > 0){
        warning("Found duplicated figure name(s) in scripts:\n",
                paste(figs_dupGen, collapse = "\n"))
    }
    if(length(figs_dupUsed) > 0){
        warning("Found plots used more than once in the same notebook::\n",
                paste(figs_dupUsed, collapse = "\n"))

    }
    figsTable$fig_exists <- file.exists(figsTable$fig)
    figsTable[order(figsTable$generated_by),]
}

check_Rmds <- function(path = ".", recursive = TRUE, full.names = TRUE){
    DT_Rmd <- data.table(
        RmdFiles = list.files(path = path, pattern = "*.Rmd",
                              recursive = recursive, full.names = TRUE) %>% 
            sort())
    DT_Rmd$html <- gsub(".Rmd$", ".html", DT_Rmd$RmdFiles)
    DT_Rmd$html_exists <- file.exists(DT_Rmd$html)
    DT_Rmd
}

get_GSEtable <- function(GSE){
    require(GEOfastq)
    require(magrittr)
    crawl_gse(GSE) %>% 
        extract_gsms() %>% 
        crawl_gsms()
}
