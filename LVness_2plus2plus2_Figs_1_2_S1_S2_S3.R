### LVness for Eberhard paper ###
### 5 vars - 3 LV not contaminated by 2 non LV ### 


rm(list = ls())




library(deSolve)

#############################################################################
### system of dif. eq. - nonLV1
### Artificial data - 	X1 and X2 LV isolated variables 
###				X3 and X4 LV exposed to every dependent variables
###				X5 and X6 non-LV variables dependent on X1 and X2

#install.packages("deSolve")
library(deSolve)

# Insert the number of time units that the simulation will run
tmax = 100
# Insert the step of the simulation
tstep = 1
t = seq(0,tmax+1,tstep) 	# time

initState_2plus2plus2=c(
	x1 = 1.2,
	x2 = .3,
	x3 = 2,
	x4 = 1,
	x5 = 1,
	x6 = 2
	) 

truePars_2plus2plus2 = c(		
   		a1 = 0.044, 
  		b11 = -0.08,
		b12 = 0.02,	
		b13 = 0,
		b14 = 0,
		b15 = 0,
		b16 = 0,
		a2 = 0.2,
		b21 = -.2,
		b22 = -.06, 
		b23 = 0,
		b24 = 0,
		b25 = 0,
		b26 = 0,
		a3 = .5,
		b31 = -.5, 
		b32 = .16,
		b33 = -.1,
		b34 = -.01,
		b35 = .3,
		b36 = -.3,
		a4 = .3,		
		b41 = -.05,
		b42 = .05,
		b43 = .2,
		b44 = -.3,
		b45 = -.4,		
		b46 = .5,		
		g51 = .156,
		f51 = -.5,
		g52 = .35,
		f52 = .6,
		g55 = -.4,
		f55 = .9,
		Vmax1 = -.14,
		k1 = .1,
		Vmax2 = .16,
		k2 = .2		
            )

Equations_2plus2plus2 <- function(t, initState, pars) 
        { 
        ### returns rate of change
        # t = model's time structure
        # initState = model initial state
        # pars = model parameters 

        with(as.list(c(initState, pars)), 
            {
            
		eq1 = x1 * ( a1 + b11 * x1 + b12 * x2 )
		eq2 = x2 * ( a2 + b21 * x1 + b22 * x2 ) 
		eq3 = x3 * ( a3 + b31 * x1 + b32 * x2 + b33 * x3 + b34 * x4 + b35 * x5 + b36 * x6 )
		eq4 = x4 * ( a4 + b41 * x1 + b42 * x2 + b43 * x3 + b44 * x4 + b45 * x5 + b46 * x6 )
		eq5 = g51 * x1^(f51) + g52 * x2^(f52) + g55 * x5^(f52)
		eq6 = Vmax1 * x1 / ( k1 + x1 ) + Vmax2 * x2 / ( k2 + x2 ) 

		dx1 = eq1
		dx2 = eq2
		dx3 = eq3
		dx4 = eq4
		dx5 = eq5
		dx6 = eq6

            return(list(c(dx1,dx2,dx3,dx4,dx5,dx6),count=c(eq1,eq2,eq3,eq4,eq5,eq6)))
            })
        }

solver <- function(times = t, initState, pars, equations = Equations) 
    {
     ### ode solves the model by integration.
     # pars = model parameters
     # equations = model equations
     # initState = model initial state
     # times = model time structure

    return(ode(times = times, y = initState, parms = pars, func = equations))
    }

### system of dif. eq. - nonLV1
#############################################################################

out_2plus2plus2 = solver(
			times = t,
			initState = initState_2plus2plus2,
			pars = truePars_2plus2plus2,
			equations = Equations_2plus2plus2)


# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\Images")     # P51
# tiff("fig222.tiff", height = 30, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(3,2))
for (i in 1:6)
	{
	plot(out_2plus2plus2[,1],out_2plus2plus2[,i+1],type='l',col="black",xlab='t',ylab=colnames(out_2plus2plus2)[i+1],cex.lab=1.5,cex.axis=1.5,lwd=3)
	}
# dev.off()

varNames = c('t',
		expression("X"[1]),
		expression("X"[2]),
		expression("X"[3]),
		expression("X"[4]),
		expression("X"[5]),
		expression("X"[6]))

parNames = c(
		expression("a"[1]),
		expression("b"[11]),
		expression("b"[12]),
		expression("b"[13]),
		expression("b"[14]),
		expression("b"[15]),
		expression("b"[16]),
		expression("a"[2]),
		expression("b"[21]),
		expression("b"[22]),
		expression("b"[23]),
		expression("b"[24]),
		expression("b"[25]),
		expression("b"[26]),
		expression("a"[3]),
		expression("b"[31]),
		expression("b"[32]),
		expression("b"[33]),
		expression("b"[34]),
		expression("b"[35]),
		expression("b"[36]),
		expression("a"[4]),
		expression("b"[41]),
		expression("b"[42]),
		expression("b"[43]),
		expression("b"[44]),
		expression("b"[45]),
		expression("b"[46]),
		expression("a"[5]),
		expression("b"[51]),
		expression("b"[52]),
		expression("b"[53]),
		expression("b"[54]),
		expression("b"[55]),
		expression("b"[56]),
		expression("a"[6]),
		expression("b"[61]),
		expression("b"[62]),
		expression("b"[63]),
		expression("b"[64]),
		expression("b"[65]),
		expression("b"[66])
		) 

axis(1, seq(-100, -50, 10), labels=labelsX)
axis(2, seq(50, 100, 10), labels=labelsY)
box()



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
##### LVness_3 - i, i+5, i+10, i+15

LVness = function (smoother_out1,suppliedSlopes=NULL)
	{
	### Function to identify if a spline data set is close to LV architecture
	# smoother_out1 # output of Smoother function
	# suppliedSlopes # sample of slopes supplied by the used. If NULL the slopes form the Smoother output will be used. Usually used to supply the true slopes

	data1 = smoother_out1$data
	DF = smoother_out1$df
	numVars = dim(data1)[2]-1
	activateWindows = FALSE
	parEst_store1 = rep(NA,numVars * (numVars+1))
	if ( numVars < 3 ) { subSampleLim = 3 } else { subSampleLim = numVars }

	STSTstartStore = rep(NA,numVars)
	depVarsSTSTStore = rep(NA,numVars)
	for (v in 2:(numVars+1))
		{
		t_length = dim(smoother_out1$splines_est)[1]
		start_low_derivs = min( max(which( (abs(smoother_out1$slopes_est[,v]) > 10^(-5)) & (abs(smoother_out1$d2X_dt2_est[,v]) > 10^(-5) ))) + 1, t_length )
		if (start_low_derivs == -Inf | start_low_derivs == t_length)
			{ 
			print('Dependent variable did not reach steady state') 
			depVarsSTSTStore[v-1] = STSTstartStore[v-1] = dim(smoother_out1$splines_est)[1]
			} else {
			print('Dependent variable is close or reach the steady state')
			#start_low_derivs_t = floor(smoother_out1$splines_est[start_low_derivs,1])
			#depVarsSTSTStore[v-1]=mean(smoother_out1$splines_est[(round(t_length*.95,0):t_length),v])
			#stst_start = min( max(which(!( (smoother_out1$splines_est[,v] < depVarsSTSTStore[v-1]*1.01) & (smoother_out1$splines_est[,v] > depVarsSTSTStore[v-1]*.99) ))) + 1, t_length )
			#stst_start_t = floor(smoother_out1$splines_est[stst_start,1])
			STSTstartStore[v-1] = start_low_derivs #max(start_low_derivs_t,stst_start_t)
			}
		}
	STST_cutoff = min(STSTstartStore,na.rm = TRUE)
	STST_cutoff_data = max(which(smoother_out1$data[,1] <= smoother_out1$splines_est[STST_cutoff,1]))
	if (STST_cutoff_data < (numVars*5+1))
		{
		print('Not enough non-STST points')
		} else {
		for (iv in 1:(min(STST_cutoff_data,dim(smoother_out1$data)[1])-(numVars*5+1)))
			{	
			if (is.null(suppliedSlopes)) 
					{
					(est1 = LV_pars_finder(
							smooth_out = smoother_out1,
							alg1_lm2 = 1, 
							data_sample_alg = iv+5*seq(0,numVars,1)		# 'random_sample'
							) )
					} else {
					smoother_out1$splines_est = smoother_out1$splines_est[smoother_out1$splines_est[,1] %in% smoother_out1$data[,1],]
					smoother_out1$slopes_est = smoother_out1$slopes_est[smoother_out1$slopes_est[,1] %in% smoother_out1$data[,1],]
					smoother_out1$d2X_dt2_est = smoother_out1$d2X_dt2_est[smoother_out1$d2X_dt2_est[,1] %in% smoother_out1$data[,1],]
					(est1 = LV_pars_finder(
							smooth_out = smoother_out1,
							supplied_slopes = suppliedSlopes,
							alg1_lm2 = 1, 
							data_sample_alg = iv+5*seq(0,numVars,1)		# 'random_sample'
							) )
					}
			parEst_store1 = cbind(parEst_store1,est1) 	
			}
		if (dim(parEst_store1)[2]!=2)
			{parEst_store1 = parEst_store1[,2:dim(parEst_store1)[2]]} 	# remove the NA line used to create matrix
		
		if (numVars < 5) 
			{par(mfrow=c(numVars,numVars+1))} else 
			{activateWindows=TRUE
			par(mfrow = c(round(sqrt(numVars+1),0)+1,round(sqrt(numVars+1),0)+1))}
			count1 = 0

		for(v in 1:dim(parEst_store1)[1])
			{
			if ( activateWindows == TRUE & count1 == numVars+1 ) {windows(); par(mfrow = c(round(sqrt(numVars+1),0)+1,round(sqrt(numVars+1),0)+1));count1 = 0}
	
			par_dist = parEst_store1[v,]
			par_slopes = par_dist[2:length(par_dist)] - par_dist[1:(length(par_dist)-1)]		
			par_slope_changes = par_slopes[2:length(par_slopes)] * par_slopes[1:(length(par_slopes)-1)]
			flips = sum(par_slope_changes<0)

			plot(	
				par_dist,
				pch=20,
				main = paste(flips,'flips'),
				ylab=rownames(parEst_store1)[v],ylim=c(-1,1)
				)
			abline(h=-1,col='lightgrey')
			abline(h=1,col='lightgrey')
			count1 = count1 + 1
			}
		return(parEst_store1)
		}
	}

##### LVness_3 - i, i+5, i+10, i+15
######################################################################################################################################

smoother_out_2plus2plus2 = Smoother(dataset = out_2plus2plus2[,1:7],data1_spline2 = 2,df = rep(100,6))	# smoothing the data

r_2plus2plus2=LVness(smoother_out1 = smoother_out_2plus2plus2)		# LV test

# Fig_S2
# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Images")     # P51
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3nd_version\\Images")     # gatech   
# tiff("pre_Fig_S2.tiff", height = 20, width = 25, units = 'cm',compression = "lzw", res = 300)
layout(matrix(c(0,rep(1,18),0,rep(1,18),0,rep(1,18),0,rep(1,18),0,rep(1,18),0,rep(1,18),
		0,rep(2,18),0,rep(2,18),0,rep(2,18),0,rep(2,18),0,rep(2,18),0,rep(2,18),
		0,rep(3,18),0,rep(3,18),0,rep(3,18),0,rep(3,18),0,rep(3,18),0,rep(3,18),
		rep(0,19),rep(0,19)), 20, 19, byrow = TRUE))			# create a frame as indicated in the matrix
layout.show(3)
par(mar = c(0.3,0,0,0),las=3)
boxplot(t(r_2plus2plus2),ylab='Estimates',xaxt='n',cex.axis=1.7)
boxplot(t(r_2plus2plus2),ylim=c(-5,5),ylab='Estimates',xaxt='n',cex.axis=1.7)
boxplot(t(r_2plus2plus2),ylim=c(-1,1),xlab='Parameters',ylab='Estimates',cex.axis=1.7,xaxt='n')
axis(side=1,at=1:42,labels=parNames,cex.axis=1.5)
# dev.off()

# tiff("pre_Fig_S2_3.tiff", height = 20, width = 25, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(1,1),mar = c(5,5,0,0))
plot(1,1,col='white',xlab = 'Parameters', ylab = 'Estimates',cex.lab=1.7)
# dev.off()


# Fig_2
# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Images")     # P51
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Images")     # gatech   
# tiff("Fig_2.tiff", height = 20, width = 40, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(2,1),mar = c(5, 5, 1, 1))
boxplot(t(r_2plus2plus2[1:21,]),ylim=c(-1,1),xlab='Parameters',ylab='Estimates',cex.lab=1.5,cex.axis=1.4,xaxt='n')
axis(side=1,at=1:21,labels=parNames[1:21],cex.axis=1.5)
boxplot(t(r_2plus2plus2[22:42,]),ylim=c(-1,1),xlab='Parameters',ylab='Estimates',cex.lab=1.5,cex.axis=1.4,xaxt='n')
axis(side=1,at=1:21,labels=parNames[22:42],cex.axis=1.5)
# dev.off()


# Fig_S1
# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Images")     # P51
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Images")     # gatech   
# tiff("Fig_S1_1.tiff", height = 30, width = 26, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(3,2),mar = c(5, 5, 1, 1),las=0)
for (i in c(1,3,5))
	{
	plot(out_2plus2plus2[,1],out_2plus2plus2[,i+1],type='l',col="black",xlab=NA,ylab=varNames[i+1],cex.lab=2,cex.axis=2,lwd=3)
	if (i == 5 | i == 6) {title(xlab="t", line=3, cex.lab=2)}	

	boxplot(t(r_2plus2plus2)[,(1:7)+((i-1)*7)],ylim=c(-1,1),cex.lab=2,cex.axis=2,xaxt='n')
	axis(side=1,at=1:7,labels=parNames[7*(i-1)+(1:7)],cex.axis=2)
	if (i != 5 & i != 6) {segments(x0=(1:7)-.3,y0=truePars_2plus2plus2[(1:7)+((i-1)*7)],x1=(1:7)+.3,col='blue',lwd=3)} else {title(xlab="Parameter", line=3, cex.lab=2)}	
	}
# dev.off()

# tiff("Fig_S1_2.tiff", height = 30, width = 26, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(3,2),mar = c(5, 5, 1, 1),las=0)
for (i in c(2,4,6))
	{
	plot(out_2plus2plus2[,1],out_2plus2plus2[,i+1],type='l',col="black",xlab=NA,ylab=varNames[i+1],cex.lab=2,cex.axis=2,lwd=3)
	if (i == 5 | i == 6) {title(xlab="t", line=3, cex.lab=2)}	

	boxplot(t(r_2plus2plus2)[,(1:7)+((i-1)*7)],ylim=c(-1,1),cex.lab=2,cex.axis=2,xaxt='n')
	axis(side=1,at=1:7,labels=parNames[7*(i-1)+(1:7)],cex.axis=2)
	if (i != 5 & i != 6) {segments(x0=(1:7)-.3,y0=truePars_2plus2plus2[(1:7)+((i-1)*7)],x1=(1:7)+.3,col='blue',lwd=3)} else {title(xlab="Parameters", line=3, cex.lab=2)}	
	}
# dev.off()



# Fig_S3
# tiff("Fig_S3.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
layout(matrix(c(0,1:7,0,8:14,0,15:21,0,22:28,0,29:35,0,36:42,0,0,0,0,0,0,0,0), 7, 8, byrow = TRUE))			# create a frame as indicated in the matrix
layout.show(42)
#par(mfrow=c(6,7),mar = c(0,0,0,0))   #,bty='n'
for (i in 1:42)
	{
	par(mar = c(0,0,0,0))
	if( i == 36 )
		{
		plot(r_2plus2plus2[i,],
			ylim=c(-1,1),
			main = NA, xlab = NA,		 			# labels
			ylab = NA,							# labels
			xaxt=NULL,yaxt=NULL					# scales
			)
		} else
		{
		plot(r_2plus2plus2[i,],
			ylim=c(-1,1),
			main = NA, xlab = NA,		 			# labels
			ylab = NA,							# labels
			xaxt='n',yaxt='n'						# scales
			)
		}
	text(x = 60,y = .9,							# position of text
		labels = parNames[i],					# the text
		col='blue',								# text color
		cex=1								# text size
		)
	}
# dev.off()








##################################################3
##### Ensemble


######################################################################################################################################
##### system of dif. eq.
##### Equations in matrix form for the pure LV system

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

smoother_out_2plus2plus2 = Smoother(dataset = out_2plus2plus2[,1:7],data1_spline2 = 2,df = rep(100,6))		# smoother

r_2plus2plus2=LVness(smoother_out1 = smoother_out_2plus2plus2)	# LV test


# Change win to check other cases
win=35

pars_est = LV_pars_finder(
				smooth_out = smoother_out_2plus2plus2,
				alg1_lm2 = 1, 
				data_sample_alg = win+c(0, 5, 10, 15, 20, 25,30) #c(1,2,5,7,9)			# c(1,6,8)	'random_sample'
				)  

initState_222 = smoother_out_2plus2plus2$splines_est[1,2:dim(smoother_out_2plus2plus2$splines_est)[2]]
names(initState_222)=colnames(out_2plus2plus2)[2:7]
cbind(pars_est)
pars_est_222_mat = Format_pars(truePars = pars_est)
out_est_222 = solveLV(
				times = seq(0,100,.1),
				initState = initState_222,
				pars = pars_est_222_mat,
				equations = Equations) 

par(mfrow=c(3,2))
for (iv in 2:7)
	{
	plot(smoother_out_2plus2plus2$data[,1],smoother_out_2plus2plus2$data[,iv],pch=9,col='darkgrey',xlab='t',ylab=varNames[iv],cex.lab=1.5,cex.axis=1.5,lwd=1)
	points(smoother_out_2plus2plus2$data[,1],smoother_out_2plus2plus2$data[,iv],pch=20,col='grey')
	points(out_est_222[,1],out_est_222[,iv],type='l',col='blue',lwd=1)
	}


#install.packages("abind")
library(abind)
fitStore = out_est_222
worked = rep(0,dim(r_2plus2plus2)[2])
initState_222 = smoother_out_2plus2plus2$splines_est[1,2:dim(smoother_out_2plus2plus2$splines_est)[2]]
names(initState_222)=colnames(out_2plus2plus2)[2:7]

for (iii in 2:dim(r_2plus2plus2)[2])
	{
	cbind(r_2plus2plus2[,iii])
	pars_est_222_mat = Format_pars(truePars = r_2plus2plus2[,iii])

	out_est_222 = NULL
	out_est_222 = try( solveLV(
					times = seq(0,100,.1),
					initState = initState_222,
					pars = pars_est_222_mat,
					equations = Equations),TRUE)

	if (class(out_est_222)!="try-error") 																					# if try does not detect an error (may create warnnings)	
		{if (dim(out_est_222)[1] == 1001)																# if out_est is the same length as the data (important to calculate errors)	
			{
			if (sum(is.nan(out_est_222))==0)																				# if there is no NAs
				{
				fitStore = abind(fitStore, out_est_222, along = 3)
				worked[iii]=1
				}
			}
		}
	}


# Fig_1
smoother_out_2plus2plus2
dim(fitStore)[3]
# setwd("C:\\Users\\dolivenca3\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Images")     # gatech
# setwd("C:\\Users\\doliv\\OneDrive\\2_2019_America\\2021\\20210427_Eberhard_slope_paper\\3rd_version\\Images")     # P51
# tiff("Fig_1.tiff", height = 20, width = 20, units = 'cm',compression = "lzw", res = 300)
par(mfrow=c(3,2),mar=c(5,5,1,1))
for (iv in 2:7)
	{
	plot(smoother_out_2plus2plus2$data[,1],smoother_out_2plus2plus2$data[,iv],pch=20,col='black',xlab='t',ylab=varNames[iv],cex.lab=1.5,cex.axis=1.5,lwd=1)
	for (v in 1:dim(fitStore)[3])
		{
		points(fitStore[,1,v],fitStore[,iv,v],type='l',col='lightblue',lwd=1)
		}
	points(smoother_out_2plus2plus2$data[,1],smoother_out_2plus2plus2$data[,iv],pch=20,col='black')
	if (iv == 2) {legend(60,1.2,legend = c('Data','ALVI-MI'),pch = c(20,NA),lty=c(NA,1),col = c('black','lightblue'),bty = "n",cex=1.3)}
	}
# dev.off()




