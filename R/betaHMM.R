#' @title The betaHMM wrapper function
#' @description A homogeneous hidden Markov model for the beta valued DNA
#'              methylation data.
#' @details A novel approach utilizing a homogeneous hidden Markov model and
#'          effectively model untransformed beta values to identify DMCs while
#'          considering the spatial correlation of the adjacent CpG sites
#' @export
#' @param data A dataframe of dimension \eqn{C \times (N \times R)} containing
#'             methylation values for \eqn{C} CpG sites from \eqn{R}
#'             treatment groups each having \eqn{N} replicates or each collected
#'             from \eqn{N} patients.
#' @param M Number of methylation states to be identified in a DNA sample.
#' @param N Number of patients or DNA sample replicates collected for each
#'          treatment group.
#' @param R Number of treatment groups (For. eg: Benign and Tumour).
#' @param seed Seed to allow for reproducibility (default = NULL).
#'
#' @return The function returns an object of the
#'         \code{\link[betaHMM:betaHMM]{betaHMM}} class
#'         which contains the following values:
#' \itemize{
#' \item function_call - The parameters passed as arguments to the function
#'                       \code{\link[betaclust:betaclust]{betaclust}}.
#' \item K - The number of hidden states identified using the BHMMs.
#' \item C - The number of CpG sites analysed using the BHMMs.
#' \item N - The number of patients or DNA replicates corresponding to each
#'           treatment group analysed using the BHMMs.
#' \item R - The number of treatment groups analysed using the BHMMs.
#' \item A - The transition matrix for the BHMM model.
#' \item tau - The initial distribution for the BHMM model.
#' \item phi - The shape parameters for the observation sequence data
#'                in the BHMM model.
#' \item log_vec - A vector containing the log-likelihood values calculated
#'                    for each iteration of the algorithm.
#' \item z - A matrix of dimension \eqn{C \times K} containing the posterior
#'              probability of each CpG site belonging to each of the
#'              \eqn{K} clusters.
#' \item hidden_states - The vector containing the estimated hidden states
#'                          for each CpG sites.
#'    }

betaHMM<-function(data,M,N,R,seed=NULL)
{
  call_function<-match.call()
  K=M^R
  C=nrow(data)
  ## Getting initialized parameters
  trained_params=initialise_parameters(data,M,N,R,seed)
  print("initialised")
  ## Baum-Welch algorithm for estimating BHMM parameters
  out_baumwelch = BaumWelch(data,trained_params,M,N,R,seed=seed)
  print("BW done")
  out_viterbi=Viterbi(data,M,N,R,out_baumwelch$tau,out_baumwelch$A,
                      out_baumwelch$phi)

  betaHMM_out<-list(function_call=call_function,K=K,C=C,N=N,R=R,
              A = out_baumwelch$A,
              tau = out_baumwelch$tau,
              phi=out_baumwelch$phi,
              log_vec=out_baumwelch$log_vec,
              z=out_baumwelch$z,
              hidden_states=out_viterbi)
  class(betaHMM_out)<-"betaHMM"
  return(betaHMM_out)
}

