# last revision: 2020/01/22


autoKrigeST = function(formula, input_data, new_data, type_stv = 'sumMetric', data_variogram = input_data, block = 0,
                     model = c("Sph", "Exp", "Gau", "Ste"),
                     kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
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

