assemble_IRF <- function(x, long = FALSE, type="IRF") {
  x_names <- names(x)
  k <- length(x_names)
  n_ahead = nrow(x[[1]])
  
  x2 <- as.data.frame(do.call("rbind", x))
  x2$impulse <- rep(x_names, each=n_ahead)
  x2$n_ahead <- rep(1:n_ahead, times = k)
  x2$type <- rep(type, times = k * n_ahead)
  x3 <- x2[, c("type",  "impulse", "n_ahead", x_names)]
  
  #long
  if(long) {
    resp_names <- colnames(x3)[4:ncol(x3)]
    x3 <- reshape(x3, 
                  direction="long", 
                  varying=resp_names,
                  times = resp_names,
                  timevar = "response",
                  v.names="value")
    x3$id <-  NULL
    
    ## the tiduverse way, so easier...
    # x4_tidy <- x3 %>%
    #   gather(response, value, -impulse, -n_ahead, -type) %>% as.tbl
  }
  
  x3
  
}

as.data.frame.varirf <- function(x, format = c("wide_resp", "long", "wide_type")) {
  format <-  match.arg(format)
  
  long <- ifelse(format %in% c("long", "wide_type"), TRUE, FALSE)
  IRF_df <- assemble_IRF(x=x$irf, type="IRF", long = long)
  if("Lower" %in% names(x)) {
    IRF_Low_df <- assemble_IRF(x$Lower, type="Lower", long = long)  
    IRF_Hig_df <- assemble_IRF(x$Upper, type="Upper", long = long)
    
    # IRF_li <- list(IRF = IRF_df, Lower = IRF_Low_df, Upper = IRF_Hig_df)
    # IRF_df <- do.call("rbind", IRF_li)
    # IRF_df_w <- Reduce(function(x, y) merge(x, y, by = c("impulse", "response", "n_ahead")), 
    #                    lapply(IRF_li, function(x) {x$type <- NULL; x}))
    # colnames(IRF_df_w) <- c("impulse", "response", "n_ahead", "IRF", "Lower", "Upper")
    
    IRF_df <- rbind(IRF_df, IRF_Low_df, IRF_Hig_df)
  }
  
  if(format == "wide_type") {
    IRF_df <- reshape(IRF_df, 
                  direction="wide", 
                  idvar = c("impulse", "n_ahead", "response"),
                  timevar = "type",
                  sep = "_") 
  }
  
  IRF_df
}

irfplot <- function(x) {
  if(inherits(x, "varirf")) df <- as.data.frame(x=x, format = "wide_type")
  irf_1_df_w2 <-  df %>% 
    mutate(impulse  = paste(impulse, "->"),
           response  = paste("->", response))

  if(all(c("value_Upper", "value_Lower") %in% colnames(irf_1_df_w2)))  {
    pl <- ggplot(aes(x = n_ahead, y=value_IRF, ymin = value_Lower, ymax = value_Upper), 
                 data= irf_1_df_w2)+
      facet_grid(impulse ~ response, switch = "y") +
      geom_smooth(stat="identity") +
      geom_hline(yintercept =0) +
      xlab("N ahead") + ylab("")
  } else {
    pl <- ggplot(aes(x = n_ahead, y=value_IRF), data= irf_1_df_w2)+
      facet_grid(impulse ~ response, switch = "y") +
      geom_line(stat="identity") +
      geom_hline(yintercept =0) +
      xlab("N ahead") + ylab("")
  }
  pl
}


if(FALSE) {
  library(tidyverse)
  library(vars)
  data(Canada)
  var.2c <- VAR(Canada, p = 2, type = "const")
  irf_1 <- irf(var.2c, runs = 1000)
  irf_1_rw <- irf(var.2c, impulse = "rw", response =  "rw", boot = TRUE, runs = 500)
  
  irf_1_df_l <- as.data.frame(x=irf_1, format = "long") %>%  as_tibble()
  irf_1_df_w1 <- as.data.frame(x=irf_1, format = "wide_resp") %>%  as_tibble()
  irf_1_df_w2 <- as.data.frame(x=irf_1, format = "wide_type") %>%  as_tibble()
  
  
  
  irfplot(irf_1)
  irfplot(irf_1_rw)
}

