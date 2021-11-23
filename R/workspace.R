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


#' Remove tmp objects
#'
#' @param collectGarbage logical. Activate garbage collector
#'
#' @return matrix By default outputs the gc() output matrix
#' @export
#'
#' @examples
rmTmp <- function(collectGarbage = TRUE){
    tmpOL <- grep(x = ls(pos = ".GlobalEnv"), pattern = "^tmp", value = TRUE)
    rm(list = tmpOL, pos = ".GlobalEnv")
    if(collectGarbage){gc()}
}

