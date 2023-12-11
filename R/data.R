#' Cross-country comparison of students' achievement
#'
#' Datasets containing average scores on math, reading, and science
#' together with standard errors for all OECD countries. These are
#' from the 2018 and 2022 editions of 
#' Program for International Student Assessment (PISA)
#' study by the Organization for Economic Cooperation and Development (OECD).
#' The average scores are over all 15-year-old students in the study.
#'
#' @format 
#' \describe{
#'   \item{jurisdiction}{country, from which data was collected}
#'   \item{math_score}{average score in math}
#'   \item{math_se}{standard error for the average score in math}
#'   \item{reading_score}{average score in reading}
#'   \item{reading_se}{standard error for the average score in reading}
#'   \item{science_score}{average score in science}
#'   \item{science_se}{standard error for the average score in science}
#' }
#' @name pisa2018
#' @source \url{https://www.oecd.org/pisa/data/}
"pisa2018"

#' @rdname pisa2018
"pisa2022"

#' Cross-country comparison of students' achievement
#'
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' New code should use `data(pisa2018)` instead.
#' 
#' Dataset containing average scores on math, reading, and science
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