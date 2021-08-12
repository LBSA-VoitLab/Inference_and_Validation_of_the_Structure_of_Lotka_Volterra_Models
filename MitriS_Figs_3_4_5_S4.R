# images for Eberhard's paper 


rm(list = ls())


# dependencies
packages <- c("deSolve","splines","gtools","fANCOVA","abind")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)



#####################################
# 	Please run LV_functions3.R	# 
#####################################

# or run the following


######################################################################################################################################
##### Smoother
#
### function to find a dataset apropriate smooth function 
# for theorical background see https://www.researchgate.net/publication/285256327_Parameter_estimation_in_canonical_biological_systems_models
#

#install.packages("splines")
library(splines)

#install.packages("fANCOVA")
library(fANCOVA)

# dataset = cbind(x1[,1],x1[,2],x2[,2]); colnames(dataset) = c('t','X1','X2')
# draw=TRUE
# data1_spline2 = 2
# smooth_type = 1 
# df = c(9,9)
# dataWeights = NULL
# splineMethod = "fmm"
# polyDegree012 = 1
# aicc1_gvc2 = 1
# span = NULL
# log_spline = FALSE

Smoother = function (
			dataset,
			draw = TRUE,
			data1_spline2 = 1, 
			smooth_type = 1,
			df = NULL,
			dataWeights = NULL, 
			splineMethod = "fmm",
			polyDegree012 = 1,
			aicc1_gvc2 = 1,
			span = NULL,
			log_spline = FALSE
			)           
	{
	### function to estimate parameters for Lotka Volterra system.
	# dataset # the data. First colunm is the timepoints, remaining colunms the dependent variables. All colunms should be named.
	# draw # TRUE draw the estimated slopes for every data point. 
	# data1_spline2 # 1 - data will be used to estimate parameters # 2 - spline samples will be used to estimate parameters 
	# smooth_type # 1 - smoothing spline that needs degrees of freedom to be defined. 2 - LOESS.
	# df # degrees of freedom for smoothing splines 2, must have 1 < df <= n. If NULL, df = n.
	# dataWeights # NULL - all datapoints have the same weights # 'heavyStart' - initial values 100 times heavier that other datapoints # it acepts handcrafted weights
	# splineMethod # method to estimate splines # "fmm" "periodic" "natural" # monotonic splines "hyman" "monoH.FC"
	# polyDegree012 # the degree of the local polynomials to be used. It can ben 0, 1 or 2.
	# aicc1_gvc2 # the criterion for automatic smoothing parameter selection: “aicc” denotes bias-corrected AIC criterion, “gcv” denotes generalized cross-validation. 
	# span # value from [0,1] to define the percentage of sample points in the moving interval
	# log_spline # TRUE uses log transform datapoints to create the splines. Splines will by in cartesian values.
	
	data_time = dataset[,1]																				# create a time var to use in splines and slopes
	numVars = dim(dataset)[2]-1																			# set the number of dependent variable	
	data_vars = matrix(rep(NA,dim(dataset)[1]*(dim(dataset)[2]-1)),ncol = numVars)										# unlisting the data
	for (iii in 1:(dim(dataset)[2]-1)) 
		{
		data_vars[,iii] = unlist(dataset[,iii+1])
		}
	colnames(data_vars) = colnames(dataset[,2:(dim(dataset)[2])])													# naming the data_vars colunms
	if (log_spline == TRUE) {data_vars = log(data_vars)}															# use log var values to calculate the splines if requested
	dataFrame = data.frame(t=dataset[,1],data_vars)							  									# create data frame for the data to use in loess  

	if ( data1_spline2 == 1 )
		{
		time_splines = data_time
		} else
		{
		time_splines = round(seq(head(data_time,1),tail(data_time,1),length.out = 10*length(data_time)),10)						# create a time sequence with a smaller resolution
		time_splines = unique(sort(c(data_time,time_splines)))
		}

	if ( !is.null(dataWeights) & length(dataWeights)==1)
		{ if ( dataWeights == 'heavyStart' ) { dataWeights = rep(1,dim(dataset)[1]);dataWeights[1]=100 } }						# set the weigths for when the initial values are 100 times heavier than the other data points

	if (is.null(df)) {df_temp = rep(NA,numVars)}

	# splines and LOESS (smoothing) 
    	splines_est = slopes_est = d2X_dt2_est = time_splines															# initiating the vars to store the results
	for (i in 1:numVars)                                                                            							# cycling thought the dependent vars
      	{
        	if (smooth_type == 1)
            	{
			if (is.null(df))
				{
				smoothSpline <- smooth.spline(data_time,data_vars[,i],w=dataWeights)
				df_temp[i] = smoothSpline$df
				} else 								# smoothing with estimated degrees of freedom 
				{smoothSpline <- smooth.spline(data_time,data_vars[,i],w=dataWeights,df=df[i])}                             		# smoothing with degrees of freedom (df)
            	smooth = predict(object = smoothSpline, x = time_splines, deriv = 0)                                 					# get the points of the fit of that linear model
            	f_of_x <- splinefun(smooth$x,smooth$y,method = splineMethod )  # "fmm" "periodic" "natural" "hyman" "monoH.FC"                                               				# create cubic spline for the linear model
            	}
	  	if (smooth_type == 2)
			{
  			loess1 = loess.as(dataFrame[,1], dataFrame[,i+1], 
						degree = polyDegree012, 
						criterion = c("aicc", "gcv")[aicc1_gvc2], 
						user.span = span, plot = F )														# create loess of the data points with span = .7
               	smooth = loess1$fit																		# store the loess fit to build the spline
                	f_of_x <- splinefun(dataFrame[,1],smooth)                                     								# create cubic spline for a dependent variable
			print(summary(loess1))
			}

		if (log_spline == FALSE)
			{		
			splines_est = cbind( splines_est, f_of_x(time_splines) )												# store spline points for all dep. variables  	
			slopes_est = cbind( slopes_est, f_of_x(time_splines, deriv = 1) ) 										# store the slopes of the data points for all dep. variables	
			d2X_dt2_est = cbind( d2X_dt2_est, f_of_x(time_splines, deriv = 2) ) 										# store the 2nd derivative estimates
			} else
			{
			splines_est = cbind( splines_est, exp(f_of_x(time_splines)) )											# store spline points for all dep. variables  	
			slopes_est = cbind( slopes_est, f_of_x(time_splines, deriv = 1) * exp(f_of_x(time_splines)) ) 						# store the slopes of the data points for all dep. variables when y is in log form. dlog(y)/dy =1/y * dy/dt <=> y * dlog(y)/dy = dy/dt	
			d2X_dt2_est = cbind( d2X_dt2_est, f_of_x(time_splines, deriv = 2) * exp(f_of_x(time_splines)) + f_of_x(time_splines, deriv = 1)^2 * exp(f_of_x(time_splines)) )  	# store the 2nd derivative estimates when y is in log form. d^2(y)/(dt)^2 = d^2(log(y))/(dt)^2 * y + (d(log(y))/dt)^2 * y 
			}
		}
	if (is.null(df)) {df = df_temp}

     	if (draw == TRUE)
            {
		par(mfrow=c(round(sqrt(numVars),0)+1,round(sqrt(numVars),0)+1)) 
		for (i in 1:numVars)                                                                            						# cycling thought the dependent vars
	      	{	
            	plot(dataset[,1],dataset[,i+1],pch=20,col="grey",ylab=colnames(data_vars)[i])                               # plot the dependent variable data

            	slopeXsize = tail(time_splines,1)*.025                                             					# find a 2.5% time interval to draw the slopes
			if (data1_spline2 == 1)			      											# draw slopes
            		{segments(x0 = time_splines - slopeXsize, x1 = time_splines + slopeXsize, y0 = dataset[,i+1] - slopeXsize * slopes_est[,i+1], y1 = dataset[,i+1] + slopeXsize * slopes_est[,i+1],col='lightgreen')} else
				{segments(x0 = time_splines - slopeXsize, x1 = time_splines + slopeXsize, y0 = splines_est[,i+1] - slopeXsize * slopes_est[,i+1], y1 = splines_est[,i+1] + slopeXsize * slopes_est[,i+1],col='lightgreen')}

            	points(splines_est[,1],splines_est[,i+1],type='l',col="darkgreen",lwd=3)                      			# plot the spline	
			points(dataset[,1],dataset[,i+1],pch=20,col="grey")

			if ( round(sqrt(numVars),0)+1 > 4 )
				{
				if (i%%16 == 0)
					{
            			windows()                                                       						# creating a new plot window, one for every dependent variable
					par(mfrow=c(4,4))
					}
				}
 			}
		}

    	return( list (
			data = dataset,
			splines_est = splines_est,			
		   	slopes_est = slopes_est,
			d2X_dt2_est = d2X_dt2_est,
			df = df
			) )
	}

##### Smoother
######################################################################################################################################



######################################################################################################################################
##### LV_par_finder
#
### function to estimate parameters for Lotka Volterra system.
# for theorical background see https://www.researchgate.net/publication/285256327_Parameter_estimation_in_canonical_biological_systems_models
#

# rm(list = ls())

# smooth_out = smoother_out
# supplied_slopes = NULL
# alg1_lm2 = 1
# data_sample_alg = 'random_sample'
# data_sample_lm = NULL
# givenParNames = NULL

LV_pars_finder = function (
				smooth_out,
				supplied_slopes = NULL,
				alg1_lm2 = 2, 
				data_sample_alg = 'random_sample',
				data_sample_lm = NULL,
				givenParNames = NULL
				)           
	{
	### function to estimate parameters for Lotka Volterra system.
	# smooth_out # please enter ypur favorite spline info
	# supplied_slopes # sample of slopes supplied by the used. If not NULL the slopes will not be calculated.
	# alg1_lm2 # 1 - the function will use a system of equations solver to find the parameters, 2 - it will use linear regretion
	# data_sample_alg # the points of the sample to be used if the solver of linear system of equations is selected. 'random_sample' - draw a random sample from the data or spline
	# data_sample_lm # the points of the sample to be used if the solver of linear system of equations is selected. When NULL it will use all sample points or spline points.
	# givenParNames # NULL - canonical par names will be asign # if the parameters have different names then the cannonical ones, you can enter then where

	dataset = smooth_out$data													# extract data from smooth_out
	data_time = dataset[,1]														# create a time var to use in splines and slopes
	numVars = dim(dataset)[2]-1													# set the number of dependent variable	
	data_vars = dataset[,2:dim(dataset)[2]]											# extract dependent vars data from smooth_out
	dataFrame = data.frame(data_time,data_vars)							  			# create data frame for the data to use in loess  
	time_splines = smooth_out$splines_est[,1]											# extracting time for splines
	splines_est = smooth_out$splines_est											# extract splines ( S(t) ) 	
	slopes_est = smooth_out$slopes_est						 						# extract slopes ( dS/dt )	
	if ( length(data_time) == length(splines_est[,1]) ) {data1_spline2 = 1} else {data1_spline2 = 2}	# see if we are using the data or spline extended time
	if (!is.null(supplied_slopes)) {slopes_est = supplied_slopes}							# if we have slopes, use it

	leftSide = vector()                                                                             	# creating the left side vector that will house the slopes / var value

    	if (alg1_lm2 == 1)
        	{
		# data sample for algebraic solution
		if (is.null(data_sample_alg)) 
			{return("Unspecified data_sample_alg")} else                                  		# if sample not defined, I define it
			{
			if (data_sample_alg[1] == 'random_sample')								# if the user wants a random sample 
				{
				data_sample_alg = sort( sample(1:length(time_splines),numVars+1) )			# take a random sample from the data or splines
				} else 
				{
				if (data1_spline2 == 2)
					{
					data_sample_alg = which(time_splines %in% data_time[data_sample_alg]) 
					}
				}
			}

		solution = vector()													# vector to store the solution
    		if (data1_spline2 == 1)
			{
			print(cbind(
					sample = data_sample_alg,
					t=time_splines[data_sample_alg],
					data_vars[data_sample_alg,]
					))
 			rightSide = cbind(rep(1,numVars+1),data_vars[data_sample_alg,])                           # build a small matrix with the values of the dependent variables values
   			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftSide = slopes_est[data_sample_alg,i+1]/data_vars[data_sample_alg,i]  		# create the left side of the system of equations with slopes / vars values
				cat('Left side ',i,"\n",leftSide,"\n")
				solution = c(solution, solve(rightSide, leftSide)) 						# solve the system of equations to find the parameter values
				}                                          		
			}
    		if (data1_spline2 == 2)
			{
			print(cbind(
					sample = data_sample_alg,
					t=time_splines[data_sample_alg],
					splines_est[data_sample_alg,2:dim(splines_est)[2]]
					))
			rightSide = cbind(rep(1,numVars+1),splines_est[data_sample_alg,2:dim(splines_est)[2]])	# build a small matrix with the values of the dependent variables values
    			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftSide = slopes_est[data_sample_alg,i+1]/splines_est[data_sample_alg,i+1]		# create the left side of the system of equations with slopes / vars values
				cat('Left side ',i,"\n",leftSide,"\n")
				solution = c(solution, solve(rightSide, leftSide)) 						# solve the system of equations to find the parameter values
				}  
			}
		cat('Right side')
		print(rightSide)
		cat("Right side determinant is ",det(rightSide),"\n")									# check the determinant of the right side of the system if equations to see if it is solvable
		cat("Right side dimensions are ",dim(rightSide)," and the rank is ",qr(rightSide)$rank,"\n","\n")	# calculate the rank of the right side of the system if equations to see if it is solvable		
        	}
    	
    	if (alg1_lm2 == 2) 
		{
		solution = vector()                                                                         	# creating vector to store the linear regretion estimates
		if (data1_spline2 == 1) 												# if we are using data
			{
			if (is.null(data_sample_lm)) {data_sample_lm = 1:length(data_time)}				# if data_sample_lm is NULL use all points in the sample
			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftside = slopes_est[data_sample_lm,i+1]/data_vars[data_sample_lm,i]			# create the leftside of the system of equations with slopes and equation var - slope/eqVar
				rightside = data_vars[data_sample_lm,] 								# create the rigthside of the system of equations with all var values 
        			solution = c(solution,lm(leftside~rightside)$coef)         					# linear regretion estimates with spline points
				}
			}
		if (data1_spline2 == 2) 												# if we are using the spline values
			{
			if (is.null(data_sample_lm)) {data_sample_lm = 1:length(time_splines)}				# if data_sample_lm is NULL use all points in the spline
			for (i in 1:numVars)                                                                      # cycling thought the dependent vars
      			{
				leftside = slopes_est[data_sample_lm,i+1]/splines_est[data_sample_lm,i+1]		# create the leftside of the system of equations with slopes and equation splin - slope/eqSpline
				rightside = splines_est[data_sample_lm,2:dim(splines_est)[2]] 				# create the rigthside of the system of equations with the spline values for all vars 
        			solution = c(solution,lm(leftside~rightside)$coef)         					# linear regretion estimates with spline points
				}
			}
		}  

	if ( !is.null(givenParNames) & length(givenParNames) == length(solution) )					# are given names for the parameters available and as long as the solutions?
		{ names(solution)=givenParNames } else										# if yes, use them
		{
		canonicalParNames = vector() 												# if no, construct cannonical pars names
		for ( i in 1:numVars )
			{
			canonicalParNames = c(canonicalParNames,paste('a',i,sep=''))					# a's
			for ( ii in 1:numVars )
				{
				canonicalParNames = c(canonicalParNames,paste('b',i,ii,sep=''))				# b's	
				}
			}
		names(solution)=canonicalParNames											# use cannonical par names
		} 
	return(solution)															# return pars_est
    	}

##### LV_par_finder
######################################################################################################################################



######################################################################################################################################
##### AlgPoinFinder

AlgPointFind = function (
				smoother_out_APF,
				dif_equations = Equations,
				matrixEq = TRUE,
				randomSearchTries = NULL,
				divideByMean = FALSE
				)
	{

	##########################################################
	### Function to find the best combination of datapoints for the albegraic method
	# smoother_out_APF # smoother results
	# dif_equations # stuff needed to run the LV solver
	# matrixEq # FALSE - parameters not in matrix form # TRUE - parameters are in matrix form
	# randomSearchTries # NULL - will check all posibilities # number - will check the given number of possibilities
	# divideByMean # FALSE - do nothing # TRUE - divide the errors by the mean of the dep var. This will balance the SSEs of the different dep. vars. if their values are very different. 

   	nVars = dim(smoother_out_APF$data)[2]-1																		# getting the number of dependent variables							

	#install.packages("gtools")
	library(gtools)
	pointComb = combinations(n=length(smoother_out_APF$data[,1]), r=nVars+1, v=1:dim(smoother_out_APF$data)[1], set=TRUE, repeats.allowed=FALSE)  	# calculating the different combinations available for the point sample
	if ( !is.null(randomSearchTries) ) { pointComb = pointComb[sample(1:dim(pointComb)[1],size = randomSearchTries),] }					# if we are doing random search, choose the random combinations to try
	worked = rep(0,dim(pointComb)[1])

	global_store = store1 = store2 = c(rep(NA,nVars+1),10^20)
	global_store_pars = rep(0,nVars * (nVars + 1))															# create the stores for the results	

	for (ii in 1:dim(pointComb)[1])																			# cycle all point combinations
		tryCatch({
			( temp1 = parEst_algPointFind = LV_pars_finder(
							smooth_out = smoother_out_APF,
							alg1_lm2 = 1, 
							data_sample_alg = pointComb[ii,]
							) )
			if (matrixEq==TRUE)
				{
				# formating the parEst_algPointFind to matrix form - use if system is in matrix form
				estPars_mat_temp = matrix(parEst_algPointFind,nrow=(-1+sqrt(1+4*length(parEst_algPointFind)))/2,byrow=TRUE)
				estPars_mat_algPointFind = list(		
									a = estPars_mat_temp[,1],
									B = estPars_mat_temp[,2:dim(estPars_mat_temp)[2]]
									)
				parEst_algPointFind = estPars_mat_algPointFind
				}

			state_algPointFind = unlist(smoother_out_APF$splines_est[1,2:(nVars+1)])
			names(state_algPointFind) = colnames(smoother_out_APF$data[,2:(nVars+1)])

			# var estimates
			out_est = NULL																								# set out_est to NULL. For previous sucessfull runs are not used when solve LV fails 
			out_est = try( solveLV(times = smoother_out_APF$data[,1], initState = state_algPointFind, pars = parEst_algPointFind, equations = dif_equations),TRUE)	# try the numerical solver	
	
			if (class(out_est)!="try-error") 																					# if try does not detect an error (may create warnnings)	
				{if (dim(out_est)[1] == length(smoother_out_APF$data[,1]))																# if out_est is the same length as the data (important to calculate errors)	
					{
					if (sum(is.nan(out_est))==0)																				# if there is no NAs
						{
						varMeans = apply(smoother_out_APF$data[,2:(nVars+1)],2,mean)													# get the means of each variable
						varMaenMatrix = matrix(rep(1,dim(smoother_out_APF$data)[1]*nVars),ncol=nVars)											# create a unitary matrix with the same dim as the data
						if ( divideByMean == TRUE ) { for (iv in 1:dim(smoother_out_APF$data)[1]) {varMaenMatrix[iv,] = varMeans} }						# if we are dividing by the means we will populate the matrix with the means of the dep. vars. matrix with the var means repeated to divide the errors so that a high value var does not dominate the errors

						error1 = sum(((smoother_out_APF$data[,2:(nVars+1)]-out_est[,2:(nVars+1)])/varMaenMatrix)^2)								# calculate the error agianst the data
						if ( error1 < store1[nVars+2] ) {store1 = c(pointComb[ii,],error1)}												# if the latest error is the smallest store it as the best 
						error2 = sum(((smoother_out_APF$splines_est[which(smoother_out_APF$splines_est[,1] %in% smoother_out_APF$data[,1]),2:(nVars+1)]-out_est[,2:(nVars+1)])/varMaenMatrix)^2)	# calculate the error agianst the splines
						if ( error2 < store2[nVars+2] ) {store2 = c(pointComb[ii,],error2) }												# if the latest error is the smallest store it as the best 
						global_store = rbind(global_store,c(pointComb[ii,],error1))														# store the best errors in the global store
						global_store_pars = cbind(global_store_pars, temp1)
						worked[ii] = 1
						}
					}
				}
		# print results and percentage of work done
		print( paste(
				paste(round(ii/dim(pointComb)[1]*100,3),'% complete || DF ='),
				paste(smoother_out_APF$df,collapse = " " ),
				paste(' || '),
				paste(round(store1,3),collapse = " " ),
 				paste(' || '),
				paste(round(store2,3),collapse = " " )
				))
		flush.console()
		})
	global_store ->> globalStore
	global_store_pars ->> globalStore_pars
	cbind(pointComb,worked) ->> work
	return( rbind(store1,store2) )
	}

##### AlgPoinFinder
######################################################################################################################################





######################################################################################################################################
##### system of dif. eq.

#install.packages("deSolve")
library(deSolve)

Format_pars = function(truePars)
	{
	truePars_temp = matrix(truePars,nrow=(-1+sqrt(1+4*length(truePars)))/2,byrow=TRUE)
	truePars_mat = list(a = truePars_temp[,1],
						B = truePars_temp[,2:dim(truePars_temp)[2]]
						)
	return(truePars_mat)
	}

Equations <- function(t, x, pars) 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(x, pars)), 
            {
		
		eq = ( a * x + x * (B %*% x) )

		dx = eq

            return(list(c(dx),count=c(eq)))
            })
        }

solveLV <- function(times = t, initState, pars, equations = Equations) 
    {
     ### ode solves the model by integration.
     # pars = model parameters
     # equations = model equations
     # initState = model initial state
     # times = model time structure

    return(ode(times = times, y = initState, parms = pars, func = equations))
    }

##### system of dif. eq.
######################################################################################################################################



##################################################################################
### data

# setwd('C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Scripts\\MitriS8_data')		# gaTech #
# setwd('C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\TestData\\Voit_data\\FW_data\\FW_data_one_looking_back')		# P51 #
# setwd('/home/dolivenca//onedrive/2_2019_America/2020/20200123_MAR/Camparison_LV_MAR/TestData/Jacob_data/Mitri_data/MitriS8')		# Linux #
# setwd('C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Scripts\\MitriS8_data')		# gaTech #

# please make the directory that contains the files your working directory
getwd()
# setwd('')

mitris_all = read.table(file = 'fta.txt', header = FALSE, sep = "", dec = ".")  

varNames = c('At','Ct','Ms','Oa')	
varColor_light = c('lightblue','orange','lightgreen','grey')
varColor = c('blue','orange3','darkgreen','black')

colnames(mitris_all) = c('Time',varNames)
mitris_all

plot(mitris_all)
cor(mitris_all[,2:5])


# replicates
mitris_rep1 = read.table(file = 'ft1.txt', header = FALSE, sep = "", dec = ".")  
mitris_rep2 = read.table(file = 'ft2.txt', header = FALSE, sep = "", dec = ".")  
mitris_rep3 = read.table(file = 'ft3.txt', header = FALSE, sep = "", dec = ".")  
mitris_reps = rbind(mitris_rep1,mitris_rep2,mitris_rep3)
mitris_reps = mitris_reps[order(mitris_reps[,1],decreasing=FALSE),]
colnames(mitris_reps) = c('Time',varNames)

par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_reps[,1],mitris_reps[,i],col=varColor[i-1],xlab='Time',ylab=varNames[i-1])
	points(mitris_all[,1],mitris_all[,i],pch=8,col=varColor[i-1],lwd=1.5)
	}

# mitris_reps_div10_6
mitris_reps_div10_6 = cbind(mitris_reps[,1], mitris_reps[,2:5]/(10^6)) 	# data
state_mitris_reps_div10_6 = mitris_reps_div10_6[1:3,2:5]			# state

### data
########################################################################




# Reduce values on y-axes by 1,000,000

###################################################
##### rescaling all divided by 10^6

mitris_all_div10_6 = cbind(mitris_all[,1], mitris_all[,2:5]/(10^6)) 	# data
state_mitris_all_div10_6 = unlist(mitris_all_div10_6[1,2:5])		# initial state

smoother_out = Smoother(
				dataset = mitris_all_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,6,5,7),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 4 4  ||  2 4 6 7 8 128300.335  ||  2 4 7 8 9 21904.32"
#  "100 % complete || DF = 4 6 5 7  ||  2 4 6 7 9 113837.948  ||  2 4 6 7 9 53476.575"

pars_est = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 4, 6, 7, 9) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris = smoother_out$splines_est[1,2:dim(smoother_out$splines_est)[2]]
names(initState_mitris)=varNames
cbind(pars_est)
pars_est_mitris_mat = Format_pars(truePars = pars_est)
out_est_mitris = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris,
				pars = pars_est_mitris_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\TestData\\Jacob_data\\Mitri_data\\300dpi_pics")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_all_div10_6[,1],mitris_all_div10_6[,i],pch=8,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	#points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_mitris[,1],out_est_mitris[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,0.07,legend = c('replicates mean','smooth','estimates'),lty = c(0,1,1),pch = c(9,NA,NA),col = c('darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black'),bty = "n",cex=1,lwd=c(NA,1,3))}
	}
# dev.off()

mitris_all_div10_6_store2 = list(splines = smoother_out$splines_est[,1:5],pars = pars_est, ests = out_est_mitris[,1:5])

##### rescaling all divided by 10^6
###################################################



########################
### mitris rep1 divided by 10^6
mitris_rep1
mitris_rep1_div10_6 = cbind(mitris_rep1[,1], mitris_rep1[,2:5]/(10^6)) 	# data
state_mitris_rep1_div10_6 = unlist(mitris_rep1_div10_6[1,2:5])		# state

smoother_out_rep1 = Smoother(
				dataset = mitris_rep1_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,5,3,4),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out_rep1,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 6 5  ||  1 2 5 7 9 165614.394  ||  1 3 5 7 9 11541.4"
#  "100 % complete || DF = 4 5 3 4  ||  1 2 5 7 9 204661.782  ||  2 4 5 7 9 1044.964"

pars_est = LV_pars_finder(
				smooth_out = smoother_out_rep1,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 4, 5, 7, 9) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris_rep1 = smoother_out_rep1$splines_est[1,2:dim(smoother_out_rep1$splines_est)[2]];
names(initState_mitris_rep1)=varNames
cbind(pars_est)
pars_est_mitris_rep1_mat = Format_pars(truePars = pars_est)
out_est1_mitris_rep1 = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris_rep1,
				pars = pars_est_mitris_rep1_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_rep1_div10_6[,1],mitris_rep1_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_div10_6[,1],mitris_all_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(smoother_out_rep1$splines_est[,1],smoother_out_rep1$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est1_mitris_rep1[,1],out_est1_mitris_rep1[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,10^9,legend = c('replicates','replicates mean','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}
# dev.off()

mitris_rep1_div10_6_store2 = list(splines = smoother_out_rep1$splines_est[,1:5],pars = pars_est, ests = out_est1_mitris_rep1[,1:5])

### mitris rep1 divided by 10^6
########################



########################
### mitris rep2 divided by 10^6
mitris_rep2
mitris_rep2_div10_6 = cbind(mitris_rep2[,1], mitris_rep2[,2:5]/(10^6)) 	# data
state_mitris_rep2_div10_6 = unlist(mitris_rep2_div10_6[1,2:5])		# state

smoother_out_rep2 = Smoother(
				dataset = mitris_rep2_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,6,4,4),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out_rep2,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 6 3  ||  2 4 7 8 9 223576.615  ||  2 4 7 8 9 5416.237"
#  "100 % complete || DF = 4 6 5 4  ||  2 4 6 8 9 126654.459  ||  2 4 6 7 9 17435.762"
#  "100 % complete || DF = 4 6 4 4  ||  2 4 7 8 9 137590.15   ||  2 4 6 8 9 9270.2"


pars_est = LV_pars_finder(
				smooth_out = smoother_out_rep2,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 4, 6, 8, 9) #c(2,4,7,8,9) 	'random_sample'
				) 

initState_mitris_rep2 = smoother_out_rep2$splines_est[1,2:dim(smoother_out_rep2$splines_est)[2]]
names(initState_mitris_rep2)=varNames
cbind(pars_est)
pars_est_mitris_rep2_mat = Format_pars(truePars = pars_est)
out_est_mitris_rep2 = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris_rep2,
				pars = pars_est_mitris_rep2_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_rep2_div10_6[,1],mitris_rep2_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_div10_6[,1],mitris_all_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(smoother_out_rep2$splines_est[,1],smoother_out_rep2$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_mitris_rep2[,1],out_est_mitris_rep2[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,10^9,legend = c('replicates','replicates mean','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}
# dev.off()

mitris_rep2_div10_6_store2 = list(splines = smoother_out_rep2$splines_est[,1:5],pars = pars_est, ests = out_est_mitris_rep2[,1:5])

### mitris rep2 divided by 10^6
########################



########################
### mitris rep3 divided by 10^6
mitris_rep3
mitris_rep3_div10_6 = cbind(mitris_rep3[,1], mitris_rep3[,2:5]/(10^6)) 	# data
state_mitris_rep3_div10_6 = unlist(mitris_rep3_div10_6[1,2:5])		# state

smoother_out_rep3 = Smoother(
				dataset = mitris_rep3_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,6,5,5),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out_rep3,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 6 5  ||  2 4 5 7 8 154350.97   ||  2 4 5 7 8 16239.488"
#  "100 % complete || DF = 4 6 5 5  ||  2 4 5 7 8 159600.389  ||  2 4 5 7 8 14137.899"

pars_est = LV_pars_finder(
				smooth_out = smoother_out_rep3,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 4, 5, 7, 8) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris_rep3 = smoother_out_rep3$splines_est[1,2:dim(smoother_out_rep3$splines_est)[2]]
names(initState_mitris_rep3)=varNames
cbind(pars_est)
pars_est_mitris_rep3_mat = Format_pars(truePars = pars_est)
out_est_mitris_rep3 = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris_rep3,
				pars = pars_est_mitris_rep3_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_rep3_div10_6[,1],mitris_rep3_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_div10_6[,1],mitris_all_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(smoother_out_rep3$splines_est[,1],smoother_out_rep3$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_mitris_rep3[,1],out_est_mitris_rep3[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,10^9,legend = c('replicates','replicates mean','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}
# dev.off()

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)

mitris_rep3_div10_6_store2 = list(splines = smoother_out_rep3$splines_est[,1:5],pars = pars_est, ests = out_est_mitris_rep3[,1:5])

### mitris rep3 divided by 10^6
########################



########################
### final plots

# Fig_S4 - splines
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\2nd_version\\Images")     # gatech
# tiff("Fig_S4.tiff", height = 24, width = 24, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(mitris_reps_div10_6[,1],mitris_reps_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_div10_6[,1],mitris_all_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(mitris_rep1_div10_6_store2[[1]][,1],mitris_rep1_div10_6_store2[[1]][,i],type='l',col=varColor_light[i-1])
	points(mitris_rep2_div10_6_store2[[1]][,1],mitris_rep2_div10_6_store2[[1]][,i],type='l',col=varColor_light[i-1])
	points(mitris_rep3_div10_6_store2[[1]][,1],mitris_rep3_div10_6_store2[[1]][,i],type='l',col=varColor_light[i-1])
	points(mitris_all_div10_6_store2[[1]][,1],mitris_all_div10_6_store2[[1]][,i],type='l',col=varColor[i-1],lwd=2)
	#if (i == 5) {legend(0,1000,legend = c('replicates','replicates mean','rep smooth','mean smooth'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor_light[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,1))}
	}
# dev.off()



# Fig_3 - fits
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\2nd_version\\Images")     # gatech
# tiff("Fig_3.tiff", height = 24, width = 24, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(mitris_reps_div10_6[,1],mitris_reps_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_div10_6[,1],mitris_all_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(mitris_rep1_div10_6_store2[[3]][,1],mitris_rep1_div10_6_store2[[3]][,i],type='l',col=varColor_light[i-1])
	points(mitris_rep2_div10_6_store2[[3]][,1],mitris_rep2_div10_6_store2[[3]][,i],type='l',col=varColor_light[i-1])
	points(mitris_rep3_div10_6_store2[[3]][,1],mitris_rep3_div10_6_store2[[3]][,i],type='l',col=varColor_light[i-1])
	points(mitris_all_div10_6_store2[[3]][,1],mitris_all_div10_6_store2[[3]][,i],type='l',col=varColor[i-1],lwd=2)
	#if (i == 5) {legend(0,1000,legend = c('replicates','replicates mean','rep smooth','mean smooth'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor_light[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,1))}
	}
# dev.off()

### final plots
########################



########################
### removing the two last points of the data

mitris_all_minus2 = mitris_all[1:(dim(mitris_all)[1]-2),]
plot(mitris_all_minus2)
cor(mitris_all_minus2[,2:5])
mitris_reps_minus2 = mitris_reps[1:(dim(mitris_reps)[1]-6),]
mitris_reps_minus2_div10_6 = cbind(mitris_reps_minus2[,1], mitris_reps_minus2[,2:5]/(10^6))


###################################################
##### rescaling all divided by 10^6

mitris_all_minus2_div10_6 = cbind(mitris_all_minus2[,1], mitris_all_minus2[,2:5]/(10^6)) 	# data
state_mitris_all_minus2_div10_6 = unlist(mitris_all_minus2_div10_6[1,2:5])		# state

smoother_out = Smoother(
				dataset = mitris_all_minus2_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,6,4,4),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 4 4  ||  1 3 4 6 7 62513.3  ||  1 3 4 6 7 21994.838"

pars_est = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(1, 3, 4, 6, 7) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris = smoother_out$splines_est[1,2:dim(smoother_out$splines_est)[2]]
names(initState_mitris)=varNames
cbind(pars_est)
pars_est_mitris_mat = Format_pars(truePars = pars_est)
out_est_mitris = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris,
				pars = pars_est_mitris_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\TestData\\Jacob_data\\Mitri_data\\300dpi_pics")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,i],pch=8,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	#points(smoother_out$splines_est[,1],smoother_out$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_mitris[,1],out_est_mitris[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,0.07,legend = c('replicates mean','smooth','estimates'),lty = c(0,1,1),pch = c(9,NA,NA),col = c('darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black'),bty = "n",cex=1,lwd=c(NA,1,3))}
	}
# dev.off()

mitris_all_minus2_div10_6_store2 = list(splines = smoother_out$splines_est[,1:5],pars = pars_est, ests = out_est_mitris[,1:5])

##### rescaling all divided by 10^6
###################################################



########################
### mitris rep1 divided by 10^6
mitris_rep1_minus2 = mitris_rep1[1:7,]
mitris_rep1_minus2_div10_6 = cbind(mitris_rep1_minus2[,1], mitris_rep1_minus2[,2:5]/(10^6)) 	# data
state_mitris_minus2_rep1_div10_6 = unlist(mitris_rep1_minus2_div10_6[1,2:5])		# state

smoother_out_rep1 = Smoother(
				dataset = mitris_rep1_minus2_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(3,4,3,5),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out_rep1,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 6 5  ||  1 2 5 7 9 165614.394  ||  1 3 5 7 9 11541.4"
#  "100 % complete || DF = 3 4 3 5  ||  1 3 4 6 7 152627.827  ||  2 3 4 6 7 3404.816"

pars_est = LV_pars_finder(
				smooth_out = smoother_out_rep1,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 3, 4, 6, 7) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris_rep1 = smoother_out_rep1$splines_est[1,2:dim(smoother_out_rep1$splines_est)[2]];
names(initState_mitris_rep1)=varNames
pars_est_mitris_rep1_mat = Format_pars(truePars = pars_est)
out_est1_mitris_rep1 = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris_rep1,
				pars = pars_est_mitris_rep1_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_rep1_minus2_div10_6[,1],mitris_rep1_minus2_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(smoother_out_rep1$splines_est[,1],smoother_out_rep1$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est1_mitris_rep1[,1],out_est1_mitris_rep1[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,10^9,legend = c('replicates','replicates mean','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}
# dev.off()

mitris_rep1_minus2_div10_6_store2 = list(splines = smoother_out_rep1$splines_est[,1:5],pars = pars_est, ests = out_est1_mitris_rep1[,1:5])

### mitris rep1 divided by 10^6
########################



########################
### mitris rep2 divided by 10^6
mitris_rep2_minus2 = mitris_rep2[1:7,]
mitris_rep2_minus2_div10_6 = cbind(mitris_rep2_minus2[,1], mitris_rep2_minus2[,2:5]/(10^6)) 	# data
state_mitris_rep2_minus2_div10_6 = unlist(mitris_rep2_minus2_div10_6[1,2:5])		# state

smoother_out_rep2 = Smoother(
				dataset = mitris_rep2_minus2_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,4,4,3),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out_rep2,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 6 3  ||  1 3 4 5 7 212822.218  ||  1 3 4 5 7 105466.004"
#  "100 % complete || DF = 4 4 4 3  ||  2 4 5 6 7 136289.07   ||  2 4 5 6 7 4480.383"

pars_est = LV_pars_finder(
				smooth_out = smoother_out_rep2,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 4, 5, 6, 7) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris_rep2 = smoother_out_rep2$splines_est[1,2:dim(smoother_out_rep2$splines_est)[2]]
names(initState_mitris_rep2)=varNames
pars_est_mitris_rep2_mat = Format_pars(truePars = pars_est)
out_est_mitris_rep2 = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris_rep2,
				pars = pars_est_mitris_rep2_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_rep2_div10_6[,1],mitris_rep2_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(smoother_out_rep2$splines_est[,1],smoother_out_rep2$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_mitris_rep2[,1],out_est_mitris_rep2[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,10^9,legend = c('replicates','replicates mean','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}
# dev.off()

mitris_rep2_minus2_div10_6_store2 = list(splines = smoother_out_rep2$splines_est[,1:5],pars = pars_est, ests = out_est_mitris_rep2[,1:5])

### mitris rep2 divided by 10^6
########################



########################
### mitris rep3 divided by 10^6
mitris_rep3_minus2 =  mitris_rep3[1:7,]
mitris_rep3_minus2_div10_6 = cbind(mitris_rep3_minus2[,1], mitris_rep3_minus2[,2:5]/(10^6)) 	# data
state_mitris_rep3_minus2_div10_6 = unlist(mitris_rep3_minus2_div10_6[1,2:5])		# state

smoother_out_rep3 = Smoother(
				dataset = mitris_rep3_minus2_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,6,6,5),
				log_spline = TRUE
				)
#AlgPointFind(
#		smoother_out_APF = smoother_out_rep3,
#		dif_equations = Equations,
#		matrixEq = TRUE
#		)
#  "100 % complete || DF = 4 6 6 5  ||  2 4 5 6 7 14867.495  ||  2 4 5 6 7 5924.04"

pars_est = LV_pars_finder(
				smooth_out = smoother_out_rep3,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 4, 5, 6, 7) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris_rep3 = smoother_out_rep3$splines_est[1,2:dim(smoother_out_rep3$splines_est)[2]]
names(initState_mitris_rep3)=varNames
pars_est_mitris_rep3_mat = Format_pars(truePars = pars_est)
out_est_mitris_rep3 = solveLV(
				times = seq(0,170,.1),
				initState = initState_mitris_rep3,
				pars = pars_est_mitris_rep3_mat,
				equations = Equations) 

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow = c(2,2))
for (i in 2:5)
	{
	plot(mitris_rep3_minus2_div10_6[,1],mitris_rep3_minus2_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(smoother_out_rep3$splines_est[,1],smoother_out_rep3$splines_est[,i],type='l',col=varColor[i-1])
	points(out_est_mitris_rep3[,1],out_est_mitris_rep3[,i],type='l',col=varColor[i-1],lwd=3)
	if (i == 5) {legend(0,10^9,legend = c('replicates','replicates mean','smooth','estimates'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,3))}
	}
# dev.off()

# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\Paper_figure")     # P51
# tiff("fig4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)

mitris_rep3_minus2_div10_6_store2 = list(splines = smoother_out_rep3$splines_est[,1:5],pars = pars_est, ests = out_est_mitris_rep3[,1:5])

### mitris rep3 divided by 10^6
########################


# images not used in the paper

# fits minus2
# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\TestData\\Jacob_data\\Mitri_data\\300dpi_pics")     # P51
# tiff("fig_fits_0_144.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(mitris_reps_minus2_div10_6[,1],mitris_reps_minus2_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(mitris_rep1_minus2_div10_6_store2[[3]][,1],mitris_rep1_minus2_div10_6_store2[[3]][,i],type='l',col=varColor_light[i-1])
	points(mitris_rep2_minus2_div10_6_store2[[3]][,1],mitris_rep2_minus2_div10_6_store2[[3]][,i],type='l',col=varColor_light[i-1])
	points(mitris_rep3_minus2_div10_6_store2[[3]][,1],mitris_rep3_minus2_div10_6_store2[[3]][,i],type='l',col=varColor_light[i-1])
	points(mitris_all_minus2_div10_6_store2[[3]][,1],mitris_all_minus2_div10_6_store2[[3]][,i],type='l',col=varColor[i-1],lwd=2)
	#if (i == 5) {legend(0,1000,legend = c('replicates','replicates mean','rep smooth','mean smooth'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor_light[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,1))}
	}
# dev.off()


# fits minus2 alternative Color schem
# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2020\\20200123_MAR\\Camparison_LV_MAR\\TestData\\Jacob_data\\Mitri_data\\300dpi_pics")     # P51
# tiff("fig_fits_0_144_altcol.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(mitris_reps_minus2_div10_6[,1],mitris_reps_minus2_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,i],pch=8,col='darkgrey',lwd=1.5)
	points(mitris_rep1_minus2_div10_6_store2[[3]][,1],mitris_rep1_minus2_div10_6_store2[[3]][,i],type='l',col='lightblue')
	points(mitris_rep2_minus2_div10_6_store2[[3]][,1],mitris_rep2_minus2_div10_6_store2[[3]][,i],type='l',col='lightblue')
	points(mitris_rep3_minus2_div10_6_store2[[3]][,1],mitris_rep3_minus2_div10_6_store2[[3]][,i],type='l',col='lightblue')
	points(mitris_all_minus2_div10_6_store2[[3]][,1],mitris_all_minus2_div10_6_store2[[3]][,i],type='l',col='blue',lwd=2)
	#if (i == 5) {legend(0,1000,legend = c('replicates','replicates mean','rep smooth','mean smooth'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor_light[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,1))}
	}
# dev.off()

### LVness
########################



########################
### compare SSEs of the 0-144 and the 0-288 until 144.

SSEs=matrix(rep(NA,8),ncol=2)
colnames(SSEs) = c('0-144','0-288')
rownames(SSEs) = c('means','rep1','rep2','rep3')

mitris_all_div10_6
mitris_all_div10_6_store2
sum((mitris_all_div10_6_store2$ests[mitris_all_div10_6_store2$ests[,1] %in% mitris_all_div10_6[,1],2:5] - mitris_all_div10_6[,2:5])^2) # SSE for all points 0 - 288
SSEs[1,2]=sum((mitris_all_div10_6_store2$ests[mitris_all_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_all_minus2_div10_6[,2:5])^2) # SSE for only 0 - 144
SSEs[2,2]=sum((mitris_rep1_div10_6_store2$ests[mitris_rep1_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_rep1_minus2_div10_6[,2:5])^2) # SSE for only 0 - 144
SSEs[3,2]=sum((mitris_rep2_div10_6_store2$ests[mitris_rep2_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_rep2_minus2_div10_6[,2:5])^2) # SSE for only 0 - 144
SSEs[4,2]=sum((mitris_rep3_div10_6_store2$ests[mitris_rep3_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_rep3_minus2_div10_6[,2:5])^2) # SSE for only 0 - 144

mitris_all_minus2_div10_6
mitris_all_minus2_div10_6_store2
SSEs[1,1]=sum((mitris_all_minus2_div10_6_store2$ests[mitris_all_minus2_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_all_minus2_div10_6[,2:5])^2)
SSEs[2,1]=sum((mitris_rep1_minus2_div10_6_store2$ests[mitris_rep1_minus2_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_rep1_minus2_div10_6[,2:5])^2)
SSEs[3,1]=sum((mitris_rep2_minus2_div10_6_store2$ests[mitris_rep2_minus2_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_rep2_minus2_div10_6[,2:5])^2)
SSEs[4,1]=sum((mitris_rep3_minus2_div10_6_store2$ests[mitris_rep3_minus2_div10_6_store2$ests[,1] %in% mitris_all_minus2_div10_6[,1],2:5] - mitris_rep3_minus2_div10_6[,2:5])^2)

SSEs

#######################



#######################
### superimpose fits from all and minus2 

# Fig_4
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\2nd_version\\Images")     # gatech
# tiff("Fig_4.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,2))
for (i in 2:5)
	{
	plot(mitris_reps_div10_6[,1],mitris_reps_div10_6[,i],pch=20,col='lightgrey',xlim=c(0,150),xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1) # all data
	points(mitris_all_div10_6[,1],mitris_all_div10_6[,i],pch=8,col='lightgrey',lwd=1.5) # all means
	points(mitris_all_div10_6_store2[[3]][,1],mitris_all_div10_6_store2[[3]][,i],type='l',col=varColor_light[i-1],lwd=2)

	points(mitris_reps_minus2_div10_6[,1],mitris_reps_minus2_div10_6[,i],pch=20,col='darkgrey',xlab='t',ylab=varNames[i-1],cex.lab=1.5,cex.axis=1.5,lwd=1) # data minus2
	points(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,i],pch=8,col='darkgrey',lwd=1.5) # means minus2
	points(mitris_all_minus2_div10_6_store2[[3]][,1],mitris_all_minus2_div10_6_store2[[3]][,i],type='l',col=varColor[i-1],lty=2,lwd=2) # mean fits minus2

	#if (i == 5) {legend(0,1000,legend = c('replicates','replicates mean','rep smooth','mean smooth'),lty = c(0,0,1,1),pch = c(20,9,NA,NA),col = c('darkgrey','darkgrey',varColor_light[i-1],varColor[i-1]),text.col=c('black','black','black'),bty = "n",cex=1,lwd=c(NA,NA,1,1))}
	}
# dev.off()



#######################################
##### mitri ensamble

### mitri ensamble 144

mitris_all_minus2_div10_6 = cbind(mitris_all_minus2[,1], mitris_all_minus2[,2:5]/(10^6)) 	# data
state_mitris_all_minus2_div10_6 = unlist(mitris_all_minus2_div10_6[1,2:5])		# state

smoother_out = Smoother(
				dataset = mitris_all_minus2_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,6,4,4),
				log_spline = TRUE
				)
AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
#  "100 % complete || DF = 4 6 4 4  ||  1 3 4 6 7 62513.3  ||  1 3 4 6 7 21994.838"

pars_est = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(1, 3, 4, 6, 7) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				) 

initState_mitris = smoother_out$splines_est[1,2:dim(smoother_out$splines_est)[2]]
names(initState_mitris)=varNames
cbind(pars_est)
pars_est_mitris_mat = Format_pars(truePars = pars_est)
out_est_mitris = solveLV(
				times = seq(0,150,.1),
				initState = initState_mitris,
				pars = pars_est_mitris_mat,
				equations = Equations) 

globalStore_pars

library(abind)
fitStore = out_est_mitris

for (iii in 2:dim(globalStore_pars)[2])
	{

	initState_mitris = smoother_out$splines_est[1,2:dim(smoother_out$splines_est)[2]]
	names(initState_mitris)=varNames
	cbind(globalStore_pars[,iii])
	pars_est_mitris_mat = Format_pars(truePars = globalStore_pars[,iii])
	out_est_mitris = solveLV(
					times = seq(0,150,.1),
					initState = initState_mitris,
					pars = pars_est_mitris_mat,
					equations = Equations) 

	fitStore = abind(fitStore, out_est_mitris, along = 3)
	}

# Fig_5
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\2nd_version\\Images")     # P51
# tiff("Fig_5.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
dim(fitStore)[3]
par(mfrow=c(2,2))
for (iv in 2:5)
	{
	plot(mitris_all_minus2_div10_6[,1],mitris_all_minus2_div10_6[,iv],pch=8,col='darkgrey',xlab='t',ylab=varNames[iv-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	for (v in 1:dim(fitStore)[3])
		{
		points(fitStore[,1,v],fitStore[,iv,v],type='l',col=varColor_light[iv-1],lwd=1)
		}
	points(fitStore[,1,1],fitStore[,iv,1],type='l',col=varColor[iv-1],lwd=3)
	}
# dev.off()



### mitri ensamble 288

mitris_all_div10_6 = cbind(mitris_all[,1], mitris_all[,2:5]/(10^6)) 	# data
state_mitris_all_div10_6 = unlist(mitris_all_div10_6[1,2:5])		# state

smoother_out = Smoother(
				dataset = mitris_all_div10_6,
				draw = TRUE,
				data1_spline2 = 2, 
				smooth_type = 1,
				df = c(4,6,5,6),
				log_spline = TRUE
				)
AlgPointFind(
		smoother_out_APF = smoother_out,
		dif_equations = Equations,
		matrixEq = TRUE
		)
#  "100 % complete || DF = 4 6 4 4  ||  2 4 6 7 8 128300.335  ||  2 4 7 8 9 21904.32"
#  "100 % complete || DF = 4 6 5 6  ||  2 4 6 7 9 101294.739  ||  2 4 6 7 9 29735.206"

pars_est = LV_pars_finder(
				smooth_out = smoother_out,
				alg1_lm2 = 1, 
				data_sample_alg = c(2, 4, 6, 7, 9) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				)  

initState_mitris = smoother_out$splines_est[1,2:dim(smoother_out$splines_est)[2]]
names(initState_mitris)=varNames
cbind(pars_est)
pars_est_mitris_mat = Format_pars(truePars = pars_est)
out_est_mitris = solveLV(
				times = seq(0,300,.1),
				initState = initState_mitris,
				pars = pars_est_mitris_mat,
				equations = Equations) 

globalStore_pars

library(abind)
fitStore = out_est_mitris

for (iii in 2:dim(globalStore_pars)[2])
	{

	initState_mitris = smoother_out$splines_est[1,2:dim(smoother_out$splines_est)[2]]
	names(initState_mitris)=varNames
	cbind(globalStore_pars[,iii])
	pars_est_mitris_mat = Format_pars(truePars = globalStore_pars[,iii])

	out_est_mitris = NULL
	out_est_mitris = try( solveLV(
					times = seq(0,300,.1),
					initState = initState_mitris,
					pars = pars_est_mitris_mat,
					equations = Equations),TRUE)

	if (class(out_est_mitris)!="try-error") 																					# if try does not detect an error (may create warnnings)	
		{if (dim(out_est_mitris)[1] == 3001)																# if out_est is the same length as the data (important to calculate errors)	
			{
			if (sum(is.nan(out_est_mitris))==0)																				# if there is no NAs
				{
				fitStore = abind(fitStore, out_est_mitris, along = 3)
				}
			}
		}

	}

# not used in paper
dim(fitStore)[3]
par(mfrow=c(2,2))
for (iv in 2:5)
	{
	plot(mitris_all_div10_6[,1],mitris_all_div10_6[,iv],pch=8,col='darkgrey',xlab='t',ylab=varNames[iv-1],cex.lab=1.5,cex.axis=1.5,lwd=1)
	for (v in 1:dim(fitStore)[3])
		{
		points(fitStore[,1,v],fitStore[,iv,v],type='l',col=varColor_light[iv-1],lwd=1)
		}
	points(fitStore[,1,1],fitStore[,iv,1],type='l',col=varColor[iv-1],lwd=3)
	}

