#' @title DMR identification from DMCs
#' @description Function to identify the DMRs for the DMCs identified using BHMM
#' @details Function to identify the DMRs for the DMCs identified using BHMM
#' @export
#' @seealso \code{\link{betaHMM}}
#' @param x A \code{\link[betaHMM:betaHMM]{betaHMM}} object.
#' @param diff_meth_cluster The clusters identified as differentially methylated
#' @param data A dataframe of dimension \eqn{C \times NR}+2 containing the
#'             chromosome, CpG site label and the methylation
#'             values for \eqn{C} CpG sites from \eqn{R} treatment groups where
#'             DNA samples are either collected from \eqn{N} patients
#'             or has \eqn{N} replicates for each treatment group and this
#'             dataset was passed as an argument to the
#'             \code{\link[betaHMM:betaHMM]{betaHMM}} function.
#' @param legacy_data A dataset containing the manifest data
#'                    from the Illumina MethylationEPIC beadchip array.
#'
#' @return A dataframe containing the following columns:
#' \itemize{
#' \item start_CpG - The starting CpG site in the particular DMR
#' \item end_CpG -  The ending CpG site in the particular DMR
#' \item DMR_size - Number of CPG sites identified in the DMR
#' \item chr_dmr - The chromosome corresponding to the CpG sites in the DMR.
#' \item map_start - MAPINFO of starting CpG site in the particular DMR
#' \item map_end -  MAPINFO of the ending CpG site in the particular DMR}
#' @importFrom stats complete.cases
dmr_identification<-function(x,diff_meth_cluster,data,legacy_data)
{

  C = x$C
  df_dmr<-cbind(data,x$hidden_states)
  block_counter <- 0
  block_start <- 0
  block_end<-0
  mat<-matrix(NA,C,3)
  block_length <- 0

  # Loop through the sequence of numbers
  for (i in 1:length(x$hidden_states)) {

    # Check if the current number is equal to 3
    if (df_dmr[i,11] == diff_meth_cluster) {

      # If this is the start of a new block, save the starting index and
      #set the block length to 1
      if (block_length == 0) {
        block_start <- i
        block_length <- 1

        # If this is part of an existing block, increment the block length
      } else {
        if(df_dmr[i,1]==df_dmr[i-1,1])
        {block_length <- block_length + 1}
        else{

          if(block_length>=2)
          {
            mat[block_counter,]<-c(block_start,i-1,block_length)
          }
          block_start <- 0
          block_length <- 0
        }


      }

      # If the block length is greater than 2,
      #increment the block counter and print the
      #starting and ending index of the block
      if (block_length >= 2) {
        block_counter <- block_counter + 1
        #cat("Block", block_counter, "starts at index", block_start,
        #"and ends at index", i, "\n")
        #mat[block_counter,]<-c(block_counter,block_start,i)
      }

      # If this is not part of a block, reset the block_start and
      #block_length variables
    } else {
      if(block_length>=2)
      {
        mat[block_counter,]<-c(block_start,i-1,block_length)
      }
      block_start <- 0
      block_length <- 0
    }
  }
  mat<-mat[stats::complete.cases(mat),]
  mat<-as.data.frame(mat)


  df_unique<-mat
  colnames(df_unique)<-c("start_CpG","end_CpG","DMR_size")
  df_unique2=df_unique
  cpg_names<-as.vector(data[,2])
  start_cpg<-sapply(df_unique2$start_CpG,function(x){cpg_names[x] })
  end_cpg<-sapply(df_unique2$end_CpG,function(x){cpg_names[x] })
  df_unique$start_CpG<-start_cpg
  df_unique$end_CpG<-end_cpg

  CHR<-sapply(df_unique$start_CpG,function(x)
    {legacy_data[legacy_data$IlmnID==x,"CHR"]})
  map_start<-sapply(df_unique$start_CpG,function(x)
    {legacy_data[legacy_data$IlmnID==x,"MAPINFO"]})
  map_end<-sapply(df_unique$end_CpG,function(x)
    {legacy_data[legacy_data$IlmnID==x,"MAPINFO"]})
  df_unique$CHR<-CHR
  df_unique$map_start<-map_start
  df_unique$map_end<-map_end

  return(df_unique)
}
