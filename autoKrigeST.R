# last revision: 2020/01/22


autoKrigeST = function(formula, input_data, new_data, type_stv = 'sumMetric', data_variogram = input_data, block = 0,
                     model = c("Sph", "Exp", "Gau", "Ste"), kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
						         fix.values = c(NA,NA,NA),
                     #remove_duplicates = TRUE,
                     newdata_mode = 'rect',
						         newdata_npoints = 1e4,
                     verbose = FALSE, GLS.model = NA,
						         tlags=0:6,
						         cutoff=2e4,
						         width=5e2,
						         aniso_method='vgm',
						         type_joint='Exp',
						         prodsum_k=0.25,
						         theoretical = FALSE,
						         start_vals = c(NA,NA,NA), miscFitOptions = list(), cores = 1, ...)
# This function performs an automatic Kriging on the data in input_data
{
  	if(inherits(formula, "STIDF") | inherits(formula, "STFDF") | inherits(formula, "STSDF"))
  	# Is someone just passes a spatialpointsdataframe, assume he/she wants to interpolate the first column with Ordinary Kriging
  	{
  		input_data = formula
  		formula = as.formula(paste(names(input_data)[1], "~ 1"))
  	}

    # Check if inpu_data and data_variogram are SpatialPointsDataFrame
    if((!inherits(input_data,"STFDF") & !inherits(data_variogram,"STFDF") & !inherits(input_data,"STIDF") & !inherits(data_variogram,"STIDF") & !inherits(input_data,"STSDF") & !inherits(data_variogram,"STSDF")))
    {
        stop(paste("\nInvalid input objects: input_data or data_variogram not of class 'ST*DF'.\n\tClass input_data: '",
                      class(input_data),"'",
                      "\n\tClass data_variogram: '",
                      class(data_variogram),"'",sep=''))
    }

  	if(as.character(formula)[3] != 1 & missing(new_data)) stop("If you want to use Universal Kriging, new_data needs to be specified \n  because the predictors are also required on the prediction spatio-temporal data elements.")
  	if("newdata" %in% names(list(...))) stop("The argument name for the prediction object is not 'newdata', but 'new_data'.")

    # TODO: remove duplicate objects in ST*DF (spatially)
    # Check if there are points or gridcells on the exact same coordinate and provide a more informative error message.
    # Points on the same spot causes the interpolation to crash.
    #if(remove_duplicates)
    #{
    #    zd = zerodist(input_data)
    #    if(length(zd) != 0)
    #    {
    #        warning("Removed ", length(zd) / 2, " duplicate observation(s) in input_data:", immediate. = TRUE)
  	#	print(input_data[c(zd), ])
    #      	input_data = input_data[-zd[, 2], ]

    #    }
    #}

    # If all the values return an informative error
    col_name = as.character(formula)[2]
    if(length(unique(input_data[[col_name]])) == 1) stop(sprintf("All data in attribute \'%s\' is identical and equal to %s\n   Can not interpolate this data", col_name, unique(input_data[[col_name]])[1]))

  	if(missing(new_data)) new_data = create_new_data.ST(input_data, form = formula, gen_mode = newdata_mode, npoints = newdata_npoints)

  	## Perform some checks on the projection systems of input_data and new_data
  	p4s_obj1 = proj4string(input_data@sp)
   	p4s_obj2 = proj4string(new_data@sp)
  	if(!all(is.na(c(p4s_obj1, p4s_obj2)))) {
  		if(is.na(p4s_obj1) & !is.na(p4s_obj2)) proj4string(input_data@sp) = proj4string(new_data@sp)
  		if(!is.na(p4s_obj1) & is.na(p4s_obj2)) proj4string(new_data@sp) = proj4string(input_data@sp)
  		#if(any(!c(is.projected(input_data@sp), is.projected(new_data@sp)))) stop(paste("Either input_data or new_data is in LongLat, please reproject.\n",
  		#									"  input_data: ", p4s_obj1, "\n",
  		#									"  new_data:   ", p4s_obj2, "\n"))
  		if(proj4string(input_data@sp) != proj4string(new_data@sp)) stop(paste("Projections of input_data and new_data do not match:\n",
  											"  input_data: ", p4s_obj1, "\n",
  											"  new_data:    ", p4s_obj2, "\n"))
  	}

    # Fit the variogram model, first check which model is used
    variogram_object = autofitVariogramST(formula = formula,
                      stf = data_variogram,
                      typestv = type_stv,
                      #autoselect.model = autoselect.model,
                      candidate_model = model,
                      tlags=tlags,
                      cutoff=cutoff,
                      width=width,
                      aniso_method=aniso_method,
                      type_joint=type_joint,
                      prodsum_k=prodsum_k,
                      theoretical = theoretical,
                      cores = cores)

    ## Perform the interpolation
    krige_result = krigeST(formula = formula,
                      data = input_data,
                      newdata = new_data,
                      modelList = variogram_object$jointSTV,
                      #block = block,
					  ...)

    krige_result@data$var.stdev <- sqrt(as.vector(krige_result@data[,'var1.pred']))

    #krige_result$var1.stdev = sqrt(krige_result$var1.var)

    # Aggregate the results into an autoKrige object
    result = list(krige_output = krige_result,var_model = variogram_object$jointSTV)
    class(result) = c("autoKrige","list")

    return(result)

}



autoKrigeST.cv <- function(stf,
                              pm = 'PM10',
                              formula = form,
                              tlags = 0:6,
                              nmax.p = 0.2,
                              LOO = FALSE,
                              tunit = NULL,
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


    kcv.vgm <- setSTI(kcv.dat, formula = formula, tlags = tlags,
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

    if (is.null(tunit)) warning('No temporal unit found. Set to the default temporal unit.')
    attr(linStAni, 'temporal unit') <- tunit

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
                  krigeST(formula = form,
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
