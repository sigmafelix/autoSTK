library(gstat)
library(spacetime)


imputeKrig.cv_dep <- function(stf, pm = 'PM10',
                              formula = as.formula(PM10~1),
                              tlags = 0:6,
                              nmax.p = 0.2, 
                              LOO = FALSE,
                              nfold = 10,
                              esti.method = 'range',
                              fit.method,
                              cutoff = 20000,
                              width = 1000,
                              variosp.type = 'Exp', 
                              variots.type = 'Ste',
                              variojo.type = 'Gau',
                              typest = 'sumMetric',
                              temp=FALSE,
                              kld.m = NULL){
  
  # target set
  stf.to <- stf
  #stf.to@data <- stf@data %>% mutate(PM10=rep(0, length(stf.to)))
  stf.to@data <- stf.to@data %>% mutate(PM10=ifelse(is.na(PM10), 0, NA_real_))
  
  ### Neither spatially nor temporally sample
  ### From here, we sample in terms of data perspective
  ### to split dataset, we use the total number of rows in target stf dataset.
  ### Shuffle indices
  index.samp = sample.int(nrow(stf.to@data), nrow(stf.to@data), replace=FALSE)
  
  #Newdata sets
  newdata.o <- stf[,,drop=FALSE]
  newdata.mis <- stf.to[,,drop=FALSE]
  
  # to subset index ... Obsolete
  nindex <- stf.to@data %>% 
    mutate(index.na = 1:length(stf.to)) %>% 
    filter(!is.na(PM10)) %>% 
    dplyr::select(index.na) %>% data.frame %>% unlist %>%  c
  # index.na <= 2
  res1 <- data.frame()
  res.name <- ifelse(LOO, paste('trial', unique(stf@sp@data[,1]), sep='.'),
                     paste('trial', formatC(1:nfold, digits = 1,flag=0), sep='.'))
  
  
  nfold <- ifelse(LOO, length(stf@sp), nfold)
  fold.k <- ceiling(length(stf@sp)/nfold)
  
  for (i in 1:nfold){
    cat("Validation group", i, "\n")
    locn <- ifelse(i==nfold, 1:(length(stf@sp) %% fold.k), 1:fold.k) + (fold.k * (i - 1)) 
    #170315 added index.pos, use it instead of locn
    index.pos = ifelse(i == nfold, 1:(length(index.samp) %% fold.k), 1:fold.k) + (fold.k * (i - 1)) 
    
    
    # This is weird. subsetting only locations instead of composite data affects the inconsistency
    # and lack of generalizability 
    # 170315
    #kcv.dat <- stf[(1:length(stf@sp))[-locn],,drop=F]
    kcv.dat = stf
    kcv.dat@data[index.pos,1] <-NA 
    
    
    kcv.vgm <- plot.sti.set(kcv.dat, formula = formula, tlags = tlags, 
                            cutoff = cutoff, width = width, wireframe = FALSE)
    kcv.vgma <- kcv.vgm
    # param settings
    ovgm.spatial = EmpVario_dist(kcv.vgma, bound = 30000, spatial=TRUE)
    ovgm.temporal = EmpVario_dist(kcv.vgma, bound = NULL, spatial=FALSE)
    
    vgm.spatial <- VarioObj_fit(orig.vgm = ovgm.spatial,
                                range = 15000, maxdist = 30000, fit.method = fit.method, plot = FALSE)
    vgm.temporal <- VarioObj_fit(orig.vgm = ovgm.temporal,
                                 range = 2, maxdist = 6, fit.method = fit.method, plot = FALSE)
    linStAni <- estiStAni(kcv.vgma, 
                          method = esti.method,
                          interval = c(100,15000),
                          spatialVgm  = vgm.spatial,
                          temporalVgm = vgm.temporal)
    
    vgm.joint <- vgm(psill = 5, variojo.type, range  = 30000, nugget = 5)
    
    attr(linStAni, 'temporal unit') <- 'hours'
    
    #######################
    print(linStAni)
    print(vgm.spatial)
    print(vgm.temporal)
    print(vgm.joint)
    print(attr(vgm.spatial, 'SSErr'))
    print(attr(vgm.temporal, 'SSErr'))
    variosp.mod = vgm.spatial
    variots.mod = vgm.temporal
    variojo.mod = vgm.joint
    
    stv.part.s = kcv.vgma %>% dplyr::filter(timelag == 0)
    stv.part.t = kcv.vgma %>% dplyr::filter(spacelag == 0)
    
    empvgm = kcv.vgma
    attr(linStAni, 'temporal unit') <- 'hours'
    
    ## sp phi, sigmasq, tausq - ts-joint - total sill and nugget,  stani
    variost.mod <- switch(typest, 
                          separable = vgmST(stModel = typest, space = variosp.mod, 
                                            time = variots.mod, sill = sill), 
                          productSum = vgmST(stModel = typest, 
                                             space = variosp.mod, time = variots.mod, k = k), 
                          productSumOld = vgmST(stModel = typest, 
                                                space = variosp.mod, time = variots.mod, 
                                                sill = psill.jo * sqrt(2), nugget = nugget), 
                          sumMetric = vgmST(stModel = typest, space = variosp.mod, time = variots.mod, 
                                            joint = variojo.mod, stAni = linStAni), 
                          simpleSumMetric = vgmST(stModel = typest,
                                                  space = variosp.mod, time = variots.mod, 
                                                  joint = variojo.mod, 
                                                  nugget = joint.nug, stAni = linStAni), 
                          metric = vgmST(stModel = typest, joint = variojo.mod, stAni = linStAni), 
                          stop(paste("model", typest, "unknown")))
    
    # strict fitModel case; Optimise the given variogram
    
    lower.stvgm = extractPar(variost.mod) * 0.5
    upper.stvgm = extractPar(variost.mod) * 1.5
    
    print(lower.stvgm);print(upper.stvgm)
    fitModel <- fit.StVariogram(kcv.vgma, 
                                variost.mod, 
                                stAni = linStAni, method = "L-BFGS-B",
                                lower = lower.stvgm,
                                upper = upper.stvgm,
                                control = list(maxit=5e3, REPORT=1)
                                )

    res1 <- rbind(res1,
                  krigeST(PM10~1, 
                          data = kcv.dat, 
                          newdata = newdata.mis,
                          fitModel, 
                          nmax = ceiling(nmax.p * nrow(stf@sp)),
                          stAni = fitModel$stAni)$var1.pred)
  }
  res1 <- data.frame(t(res1[colSums(!is.na(res1)) > 0]))
  res <- list(res1)
  return(res)
}
