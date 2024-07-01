"""
Pricing problem
"""
function Pricing(modelp,Thetamax,CoverageSatLat,CoverageSatLong,dualvar)

# Variables de linéarisation
	@variable(modelp, Xi[1:1, 1:NUM_TIME, 1:NUM_PIXEL], Bin)
# Fonction objectif : Minimize the number of active satellite
	@objective(modelp, Min, sum(dualvar[j,k] * (sum(Xi[i,p,j] for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k) for i in [1])) for j in 1:NUM_PIXEL for k=1:n_0))
# Somma su indici delle orbite
	@constraint(modelp, sum((Xi[i,p,j]) for i=[1] for p in 1:NUM_TIME for j in 1:NUM_PIXEL) >= 1)
	@constraint(modelp, [i=[1],p=1:NUM_TIME,j=1:NUM_PIXEL], Thetamax[i,p,j] >= CoverageSatLat[i,p,j] - pi*(1-Xi[i,p,j]))
	@constraint(modelp, [i=[1],p=1:NUM_TIME,j=1:NUM_PIXEL], Thetamax[i,p,j] >= CoverageSatLong[i,p,j] - 2*pi*(1-Xi[i,p,j]))
	optimize!(modelp)
	return modelp
end

"""
Restrained Master Problem
"""
function RMS(model,INDEX,Thetamax,CoverageSatLat,CoverageSatLong)
	sizeIndex = length(INDEX)
# Binary activation of satellite
	@variable(model, 0 <= activSat[1:sizeIndex] <= 1)
	@variable(model, 0 <= Xi[1:sizeIndex, 1:NUM_TIME, 1:NUM_PIXEL] <= 1)
# Fonction objectif : Minimize the number of active satellite
	@objective(model, Min, sum(activSat[i] for i in 1:sizeIndex))
	@constraint(model, observable[j=1:NUM_PIXEL,k=1:n_0], -sum(Xi[i,p,j] for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k) for i in 1:sizeIndex) <= -1)
	@constraint(model, [i=1:sizeIndex,p=1:NUM_TIME,j=1:NUM_PIXEL], Thetamax[i,p,j] >= CoverageSatLat[i,p,j] - pi*(1-Xi[i,p,j]))
	@constraint(model, [i=1:sizeIndex,p=1:NUM_TIME,j=1:NUM_PIXEL], Thetamax[i,p,j] >= CoverageSatLong[i,p,j] - 2*pi*(1-Xi[i,p,j]))
	@constraint(model, [i=1:sizeIndex, p=1:NUM_TIME,j=1:NUM_PIXEL], Xi[i,p,j] <= activSat[i])

	optimize!(model)
	status = termination_status(model)
	#println(status)
	#satellites = value.(model[:activSat])
	#Xis = value.(model[:Xi])

	#println("NUM_SAT : ", sum(satellites[i] for i in 1:sizeIndex))
	#println("NUM_OBS : " , sum((1-Xis[i,p,j]) for i in 1:sizeIndex for p in 1:NUM_TIME for j in 1:NUM_PIXEL))
#=
# Convert INDEX vector in index of Keplerian's parameters
	IndexAltitude,Indexinclination,IndexRaan,IndexMeanAnom = ExtractKepParamFromINDEX(INDEX)

# Show a binary array for each target and time step to see if solution is valide
	TimeSeen_track = Validation(altitudeSet,InclinaisonSet,RAANSet,MeanAnomalySet,IndexMeanAnom,IndexRaan,Indexinclination,IndexAltitude,inst)

	#global observability = fill(0.0, 1:NUM_PIXEL, 1:n_0)
	observability = 0.0
	flag_observed = 1.0

	#println("observed period(s)")
	for j in 1:NUM_PIXEL
		for k in 1:n_0
			observability = 0.0
			for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k)
				observability += TimeSeen_track[j,p]
			end
			if observability < 1
				flag_observed = 0.0
				break
			else
				println(k)
			end
		end
	end
	#println("")
=#
	dualvar = dual.(observable)
	return model,dualvar
end

"""
Original SCDP to solve to verify the solution
"""
function Original(model1,INDEX,Thetamax,CoverageSatLat,CoverageSatLong)

	sizeIndex = length(INDEX)
# Binary activation of satellite
	@variable(model1, activSat[1:sizeIndex], Bin)
# Variables de linéarisation
	@variable(model1, Xi[1:sizeIndex, 1:NUM_TIME, 1:NUM_PIXEL], Bin)
# Fonction objectif : Minimize the number of active satellite
	@objective(model1, Min, sum(activSat[i] for i in 1:sizeIndex))
# Somma su indici delle orbite
	@constraint(model1, [j=1:NUM_PIXEL,k=1:n_0], sum((Xi[i,p,j]) for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k) for i in 1:sizeIndex) >= 1)
	@constraint(model1, [i=1:sizeIndex], sum((Xi[i,p,j]) for p in 1:NUM_TIME for j in 1:NUM_PIXEL) - 0.5 <= NUM_TIME*NUM_PIXEL*activSat[i])
	@constraint(model1, [i=1:sizeIndex], 0.5 - sum((Xi[i,p,j]) for p in 1:NUM_TIME for j in 1:NUM_PIXEL) <= NUM_TIME*NUM_PIXEL*(1-activSat[i]))
	##@constraint(model1, coplat[i=1:sizeIndex,p=1:NUM_TIME,j=1:NUM_PIXEL], Xi[i,p,j] => {Thetamax[i,p,j] >= CoverageSatLat[i,p,j]})
	##@constraint(model1, coplong[i=1:sizeIndex,p=1:NUM_TIME,j=1:NUM_PIXEL], Xi[i,p,j] => {Thetamax[i,p,j] >= CoverageSatLong[i,p,j]})
	@constraint(model1, [i=1:sizeIndex,p=1:NUM_TIME,j=1:NUM_PIXEL], Thetamax[i,p,j] >= CoverageSatLat[i,p,j] - pi*(1-Xi[i,p,j]))
	@constraint(model1, [i=1:sizeIndex,p=1:NUM_TIME,j=1:NUM_PIXEL], Thetamax[i,p,j] >= CoverageSatLong[i,p,j] - 2*pi*(1-Xi[i,p,j]))

	optimize!(model1)
	status = termination_status(model1)
	#println(status)

	obs = zeros(size(INDEX)[1])
	INDEX_FINAL = []

	indx = 0
	XiTempModel1 = value.(model1[:Xi])
	for i in 1:sizeIndex
		indx += 1
		obs[indx] = sum(XiTempModel1[i,p,j] for p in 1:NUM_TIME for j in 1:NUM_PIXEL)
		if obs[indx] >= 1
			append!(INDEX_FINAL,INDEX[i])
		end
	end
	#println("NUM_OBS : " , obs)
	#println("NUM_SAT : ", size(INDEX_FINAL)[1])

	return model1,INDEX_FINAL,sum(obs),XiTempModel1
end
