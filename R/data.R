#' Cross-country comparison of students' achievement
#'
#' A dataset containing average scores on math, reading, and science
#' together with standard errors for all OECD countries. These are
#' from the 2018 Program for International Student Assessment (PISA)
#' study by the Organization for Economic Cooperation and Development (OECD).
#' The average scores are over all 15-year-old students in the study.
#'
#' @format A data frame with 37 rows and 6 variables:
#' \describe{
#'   \item{jurisdiction}{country, from which data was collected}
#'   \item{math_score}{average score in math}
#'   \item{math_se}{standard error for the average score in math}
#'   \item{reading_score}{average score in reading}
#'   \item{reading_se}{standard error for the average score in reading}
#'   \item{science_score}{average score in science}
#'   \item{science_se}{standard error for the average score in science}
#' }
#' @source \url{https://www.oecd.org/pisa/data/}
"pisa"

#' Income of parents and children
#'
#' An artificial dataset containing income of children and their parents together
#' with some information about them. 
#'
#' @format A data frame with 3894 rows and 4 variables:
#' \describe{
#'   \item{c_faminc}{Family income of a child}
#'   \item{p_faminc}{Family income of parent}
#'   \item{gender}{Gender}
#'   \item{race}{Race: hisp (Hispanic), black or neither}
#' }
"parent_child_income"