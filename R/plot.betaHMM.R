globalVariables(c("Patient_Sample","label","Uncertainty","M","beta_value"))
#' @title Plots for visualizing the betaHMM class object
#' @description Visualise a \code{\link[betaHMM:betaHMM]{betaHMM}} clustering solution by plotting the fitted and kernel density estimates and the uncertainty.
#' @details The fitted density estimates can be visualized under the optimal clustering solution by specifying what = "fitted density" and kernel density estimates under the optimal clustering solution by specifying what = "kernel density".
#' @rdname plot.betaHMM
#' @exportS3Method plot betaHMM
#' @seealso \code{\link{betaHMM}}
#' @param x A \code{\link[betaHMM:betaHMM]{betaHMM}} object.
#' @param what The different plots that can be obtained are either
#'             "fitted density","kernel density" or
#'             "uncertainty" (default = "fitted density").
#' @param data A dataframe of dimension \eqn{C \times NR} containing methylation
#'             values for \eqn{C} CpG sites from \eqn{R} treatment groups where
#'             DNA samples are either collected from \eqn{N} patients
#'             or has \eqn{N} replicates for each treatment group and this
#'             dataset was passed as an argument to the
#'             \code{\link[betaHMM:betaHMM]{betaHMM}} function.
#' The data is not required as an input when generating "uncertainty" plot and
#'  the default has been set as "NULL". The data needs to be passed as an
#'  argument to this function when generating either
#'  "fitted density" or "kernel density" plots.
#' @param treatment_group The names of the different treatment groups.
#'  If no value is passed then default values of sample names, e.g. Sample 1,
#'  Sample 2, etc are used as legend text (default = NULL).
#' @param title The title that the user wants to display.
#' If no title is to be displayed the default is "NULL".
#' @param ... Other graphics parameters.
#'
#' @return This function displays the following plots as requested by the user:
#' \itemize{
#' \item fitted density estimates - Plot showing the fitted density estimates of the clustering solution under the optimal model selected.
#' \item kernel density estimates - Plot showing the kernel density estimates of the clustering solution under the optimal model selected.
#' \item uncertainty -  A boxplot showing the uncertainties in the optimal clustering solution.
#' }
#'
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom stats C
#' @importFrom scales seq_gradient_pal
#'
plot.betaHMM<-function(x,
                       what = "fitted density",
                       data = NULL,
                       treatment_group = NULL,
                       title = NULL,...)
{
  #UseMethod("plot")
  object <- x

if(is.null(title))
{
  txt=""
}else{
  txt=title
}
if(is.null(treatment_group))
{
  R=object$R
  treatment_group=sapply(1:R, function(x) paste0("Sample ",x))
}
if(what == "kernel density")
{
  if(is.null(data)){
    plot_graph=NULL
    warning("data argument cannot be NULL for generation of kernel density plot.", call. = FALSE)
  }else
  {


    # if(object$optimal_model == "K.." || object$optimal_model == "KN.")
    # {
    #   #data<-object$optimal_model_results$data
    #   column_len=ncol(data)
    #   if(column_len>object$N)
    #   {
    #     call_data=data[,1:object$N]
    #   }else{
    #     call_data=data
    #   }
    #   data_ggplot<-as.data.frame(call_data)
    #   #data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
    #   data_ggplot$mem_final<-as.factor(object$optimal_model_results$classification)
    #   colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"
    #   plot_graph<-ggplot2::ggplot(data_ggplot,ggplot2::aes(x=data_ggplot[,pn], fill=Cluster))+
    #     ggplot2::geom_density(alpha=0.6)+
    #     ggplot2::labs(x="Beta value", y="Density",
    #                   #title=paste0("Density estimates for ",object$optimal_model,"  clustering solution"),
    #                   title=txt,
    #                   fill ="Cluster")
    #   if(object$K==3)
    #   {
    #     colours<-c("chartreuse3","magenta","cyan3")
    #     plot_graph<-plot_graph+
    #       ggplot2::scale_fill_manual(values=colours)
    #   }
    #
    #   p.data <- ggplot2::ggplot_build(plot_graph)$data[[1]]
    #
    #
    #   p.text <- lapply(split(p.data, f = p.data$group), function(df){
    #     df[which.max(df$scaled), ]
    #   })
    #   p.text <- do.call(rbind, p.text)
    #   p.text$prop=p.text$n/(sum(p.text$n))
    #
    #
    #   plot_graph<-plot_graph + ggplot2::annotate('text', x = p.text$x,
    #                                              y = p.text$y,
    #                                              label = sprintf('%.3f',
    #                                                              p.text$prop),
    #                                              colour=p.text$fill,
    #                                              vjust = 0)
    #
    #
    # }else
    # {
      column_len=ncol(data)
      if(column_len==(object$N*object$R))
      {
        call_data=data
      }else if(column_len>(object$N*object$R)){
        call_data=data[,1:(object$N*object$R)]
      }else{
        call_data=NULL
      }
      #data<-object$optimal_model_results$data
      if(is.null(call_data))
      {
        plot_graph=NULL
        warning("K.R considers only multiple samples. Please pass DNA
                methylation values for more than 1 sample.", call. = FALSE)
      }
      else{
        data_ggplot<-as.data.frame(call_data)
        #data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
        data_ggplot$mem_final<-as.factor(object$hidden_states)
        colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"
        cols=ncol(data_ggplot)
        rows=nrow(data_ggplot)
        data_matrix<-as.matrix(data_ggplot[,1:(cols-1)])
        data_new<-as.vector(data_matrix)
        col_names<-colnames(data_ggplot)
        col_len<-length(col_names)
        Cluster<-vector()
        Patient_sample=vector()
        for(i in 1:(col_len-1))
        {
          temp=gsub("_"," ",col_names[i])
          ps_names<-rep(temp,rows)
          Patient_sample<-c(Patient_sample,ps_names)
          Cluster<-c(Cluster,data_ggplot[,cols])
        }
        data_plot<-cbind(data_new,Cluster,Patient_sample)
        data_plot<-as.data.frame(data_plot)
        colnames(data_plot)<-c("beta_value","Cluster","Patient_Sample")
        data_plot$beta_value<-as.numeric(data_plot$beta_value)
        color_length<-col_len-1
        colours<-scales::seq_gradient_pal(low="#FFC20A",
                                          high="#0C7BDC",
                                          space = "Lab"
                                          )(1:color_length/color_length)

        plot_graph<-ggplot2::ggplot(data_plot)+
          ggplot2::geom_density(aes(x=beta_value,color=Patient_Sample))+
          ggplot2::xlab("Beta Value")+
          ggplot2::ylab("Density")+
          ggplot2::scale_color_manual("Treatment group",values=colours)+
          ggplot2::facet_wrap(~Cluster,scales = "free_y"
          )+
          ggplot2::theme(axis.title.x = ggplot2::element_text(size=10),
                         axis.title.y = ggplot2::element_text(size=10)) +
          ggplot2::ggtitle(txt)
        #ggplot2::ggtitle("Density estimates for K.R clustering solution")

        cluster_size=table(object$hidden_states)
        f_labels<-data.frame(Cluster=seq(1,length(cluster_size),by=1),
                             label=as.vector(round((cluster_size/object$C),3)))
        plot_graph<-plot_graph+
          ggplot2::geom_text(x = 0.2, y = 1, ggplot2::aes(label = label),
                             data = f_labels)

      }
    }

  # }


}
if(what=="fitted density")
{
  # if(is.null(data)){
  #   plot_graph=NULL
  #   warning("data argument cannot be NULL for generation of fitted density
  #           plot.", call. = FALSE)
  # }else
  # {
    # if(object$optimal_model == "K.." || object$optimal_model == "KN.")
    # {
    #
    #   #data_x=sort(object$optimal_model_results$data[,pn])
    #   data_x=sort(data[,pn])
    #   K=object$K
    #   prop=object$optimal_model_results$tau
    #   data_th_plot<-matrix(data=NA,nrow=1,ncol=3)
    #   data_th_plot<-as.data.frame(data_th_plot)
    #   colnames(data_th_plot)<-c("beta_value","density","Cluster")
    #   alpha=object$optimal_model_results$alpha
    #   delta=object$optimal_model_results$delta
    #   if(object$optimal_model == "K..")
    #   {
    #     for(i in 1:K)
    #     {
    #       beta_value=data_x
    #       density=prop[i]*stats::dbeta(data_x,alpha[i],delta[i])
    #       Cluster<-rep(i,length(data_x))
    #       temp<-cbind(beta_value,density,Cluster)
    #       data_th_plot<-rbind(data_th_plot,temp)
    #     }
    #   }else if(object$optimal_model == "KN.")
    #   {
    #     for(i in 1:K)
    #     {
    #       beta_value=data_x
    #       density=prop[i]*stats::dbeta(data_x,alpha[i,pn],delta[i,pn])
    #       Cluster<-rep(i,length(data_x))
    #       temp<-cbind(beta_value,density,Cluster)
    #       data_th_plot<-rbind(data_th_plot,temp)
    #     }
    #   }
    #   data_th_plot<-as.data.frame(data_th_plot)
    #   data_th_plot<-data_th_plot[-1,]
    #   data_th_plot$Cluster<-as.factor(data_th_plot$Cluster)
    #   #txt=""
    #   #colours<-c("chartreuse3","magenta","cyan3")
    #   plot_graph<-ggplot2::ggplot(data_th_plot)+
    #     ggplot2::geom_line(ggplot2::aes(beta_value,density,color=Cluster),
    #                        linetype = "solid")+
    #     ggplot2::labs(x="Beta value",y="Density",title=txt,
    #                   color ="Cluster")
    #   # +
    #   #   ggplot2::scale_color_manual(values=colours)
    #   if(K==3)
    #   {
    #     colours<-c("chartreuse3","magenta","cyan3")
    #     plot_graph<-plot_graph+
    #       ggplot2::scale_color_manual(values=colours)
    #   }
    #   p.data <- ggplot2::ggplot_build(plot_graph)$data[[1]]
    #
    #
    #   p.text <- lapply(split(p.data, f = p.data$group), function(df){
    #     df[which.max(df$y), ]
    #   })
    #   p.text <- do.call(rbind, p.text)
    #   p.text$prop=prop
    #   plot_graph<-plot_graph + ggplot2::annotate('text', x = p.text$x,
    #                                              #y = p.text$y,
    #                                              y=0.2,
    #                                              label = sprintf('%.3f',
    #                                                              p.text$prop),
    #                                              colour=p.text$colour,
    #                                              vjust = 0)
    # }
    # else{
      vec_C=1001
      R<-object$R
      K<-object$K
      N<-object$N
      C<-object$C
      alpha<-t(object$phi$sp_1)
      delta<-t(object$phi$sp_2)
      tau<-round(as.vector(table(object$hidden_states)/
                           length(object$hidden_states)),3)

      density_vec<-vector()
      cluster_vec<-vector()
      sample_vec<-vector()
      beta_vec<-vector()
      vec_x<-seq(0.001, 0.999, length=vec_C)
      #treatment_group<-c("Sample A","Sample B")

      for(i in 1:R)
      {
        for(j in 1:K)
        {
          #j=1
          tmp_vec<-sapply(vec_x,function(x) {tau[j]*
                                            (stats::dbeta(x,alpha[j,i]
                                            ,delta[j,i]))})
          beta_vec<-c(beta_vec,vec_x)
          density_vec<-c(density_vec,tmp_vec)
          cluster_vec<-c(cluster_vec,rep(j,times=length(tmp_vec)))
          sample_vec<-c(sample_vec,rep(treatment_group[i],times=length(tmp_vec)))

        }
      }

      df_new_tmp<-as.data.frame(cbind(beta_vec,density_vec,cluster_vec,
                                      sample_vec))
      df_new_tmp$sample_vec<-as.factor(df_new_tmp$sample_vec)
      df_new_tmp$cluster_vec<-as.factor(df_new_tmp$cluster_vec)
      df_new_tmp$beta_vec<-as.numeric(df_new_tmp$beta_vec)
      df_new_tmp$density_vec<-as.numeric(df_new_tmp$density_vec)
      color_length<-R
      colours<-scales::seq_gradient_pal(low="#FFC20A",
                                        high="#0C7BDC",space =
                                          "Lab")(1:color_length/color_length)
      cluster_size=table(object$hidden_states)
      plot_graph<-ggplot2::ggplot(df_new_tmp,ggplot2::aes(x=beta_vec,
                                                          y=density_vec,
                                                          color=sample_vec))+
        ggplot2::geom_line()+
        ggplot2::scale_color_manual(values=colours)+
        ggplot2::facet_wrap(~cluster_vec,scales = "free_y"
        )+ ggplot2::labs(color="DNA sample", x="Beta value", y="Density")+
        ggplot2::ggtitle(txt)
      f_labels<-data.frame(cluster_vec=as.factor(seq(1,K,by=1)),
                           label=as.vector(round((cluster_size/object$C),3)))
      plot_graph<-plot_graph+
        ggplot2::geom_text(data = f_labels, ggplot2::aes(x = 0.2, y = 0.1,
                                                         label = label,
                                                         color=NA),
                                      show.legend = F,fontface="bold" )

    }

  # }

# }
# if(threshold==TRUE)
# {
#   if(!is.null(plot_graph)){
#     if(object$optimal_model == "K..")
#     {
#       th_plot<-object$optimal_model_results$thresholds
#       ano_th<-object$optimal_model_results$thresholds-.02
#     }else{
#       th_plot<-object$optimal_model_results$thresholds$threholds[,pn]
#       ano_th<-object$optimal_model_results$thresholds$threholds[,pn]-.02
#     }
#     num=sort(p.text$y)
#     ano_y=num[length(num)]-0.1
#     plot_graph<-plot_graph+ggplot2::geom_vline(xintercept = th_plot,linetype="dotted")+ggplot2::annotate("text",x=ano_th,y=ano_y,label=th_plot,angle=90)
#   }
# }

if(what == "uncertainty")
{
  labels<-c(max_uc="Maximum uncertainty")
  classification_final=apply(object$z,1,max)
  uncertainty=1-classification_final
  tau=(table(object$hidden_states))/object$C
  unc_df<-cbind(uncertainty,object$hidden_states)
  unc_df<-as.data.frame(unc_df)
  colnames(unc_df)<-c("Uncertainty","Cluster")
  unc_df$Cluster<-as.factor(unc_df$Cluster)
  unc_df_sorted<-unc_df[order(unc_df$Cluster),]
  max_unc=1-1/(M^R)
  h=max_unc+0.015
  max_uncertainty <- data.frame(yintercept=h, max_uncertainty=factor(h))
  plot_graph<-ggplot2::ggplot(unc_df_sorted, ggplot2::aes(x=Cluster,
                                                          y=Uncertainty,
                                                          color=Cluster)) +
    ggplot2::geom_boxplot()+
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank()
    )+ggplot2::labs(color="Cluster number")+
    ggplot2::xlab("Cluster number")+
    #ggplot2::theme(legend.position = "none")+
    ggplot2::ggtitle(txt)+
    #ggplot2::ggtitle("Boxplot for uncertainties in clustering solution")+
    ggplot2::coord_cartesian(ylim = c(0, 1))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=max_unc,
                                     linetype="Maximum uncertainty"),
                        color="black")+
    ggplot2::scale_linetype_manual(name="",
                                   values = 2,guide =
                                     ggplot2::guide_legend(override.aes =
                                                             list(color =
                                                                    "black")))



}

  plot_graph

# if(!is.null(plot_graph))
# {
#   if(plot_type=="ggplot")
#     plot_graph
#   else
#     plotly::ggplotly(plot_graph)
# }
}

