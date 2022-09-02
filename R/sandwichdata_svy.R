#' Example Dataset: sandwichdata_svy
#'
#' This sample data set was simulated to create variables
#' that are similar to ones we might see in a real data
#' example from a complex survey design where we want to apply
#' the 2-stage approach of regression calibration.
#' We use this data set to generate the workflow in the
#' main package vignette for the case of a complex survey.
#'
#' See package vignettes for more details.
#'
#' @docType data
#'
#' @name sandwichdata_svy
#' @format \code{sandwichdata_svy:} a data frame with 10000 rows and 9 columns
#' \describe{
#'  \item{\code{ID}}{a unique ID variable for each subject in the data frame}
#'  \item{\code{PSUid}}{the primary sampling unit (PSU) or cluster ID variable}
#'  \item{\code{strat}}{the variable specifying the strata of the survey design}
#'  \item{\code{myweights}}{the sampling weights reflecting unequal probability of selection into the sample}
#'  \item{\code{xstar}}{the error-prone, continuous exposure variable, available on all cohort study subjects (e.g. self-reported energy intake)}
#'  \item{\code{xstarstar}}{the exposure variable prone to classical measurement error, available only on a subset (e.g. recovery biomarker for energy intake)}
#'  \item{\code{v}}{an indicator variable for whether subjecs are in the subset used to fit the stage 1 model}
#'  \item{\code{z}}{a precisely-measured continuous covariate in the data set, also assumed to be related to the outcomes of interest}
#'  \item{\code{y}}{a binary outcome variable in the data set (e.g. hypertension) to be used in a logistic regression model}
#'  \item{\code{Time}}{failure time in years, to be used in a Cox model}
#'  \item{\code{delta}}{an indicator of whether failure or right-censoring occured at the end of the time period, to be used in a Cox model}
#'  }
NULL
