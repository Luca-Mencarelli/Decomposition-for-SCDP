push!(LOAD_PATH, "$(pwd())/SCDP")
using JuMP, SCDP, StaticArrays, SatelliteToolbox, Random, JLD2, CPLEX, DataStructures, Ipopt, Distributions
include("SCDP/SCDP/src/constantes_physiques.jl")
include("RMPv3a.jl")

Random.seed!(0)

function calcul_angle_ECI_ECEF(date::Real)
	mat_ECI_ECEF = r_eci_to_ecef(TOD(), PEF(), date)
	return atan(-mat_ECI_ECEF[2, 1], mat_ECI_ECEF[1, 1])
 end

global RAYON = 6378136.3 # m
global angle0 = calcul_angle_ECI_ECEF(time_zero_simulation)

function InitialSolution(modelf,Latitude,Longitude,time)
	global NUM_SATELLITE = 1
	global NUM_PIX = 1
	global cputime = 0

	register(modelf, :%, 2, %; autodiff = true)

	# Variables de linéarisation
	@variable(modelf, Thetamax[1:NUM_SATELLITE])

	# Altitude de l'orbite [Km]
	@variable(modelf, altitude[1:NUM_SATELLITE] >= RAYON)
	@variable(modelf, 0 <= activation_altitude[1:NUM_SATELLITE,1:length(altitudeSet)] <= 1)

	@constraint(modelf, [i=1:NUM_SATELLITE,a=1:length(altitudeSet)], activation_altitude[i,a]*(1-activation_altitude[i,a]) <= 10e-5)
	@constraint(modelf, [i=1:NUM_SATELLITE,a=1:length(altitudeSet)], activation_altitude[i,a]*(1-activation_altitude[i,a]) >= -10e-5)

	# inclinaison de l'orbite
	@variable(modelf, 0 <= sin_inclinaison[1:NUM_SATELLITE] <= 1)
	@variable(modelf, -1 <= cos_inclinaison[1:NUM_SATELLITE] <= 1)

	# Noeud ascendant
	@variable(modelf, -1 <= sin_noeudAscendant[1:NUM_SATELLITE] <= 1)
	@variable(modelf, -1 <= cos_noeudAscendant[1:NUM_SATELLITE] <= 1)

	# Anomalie moyenne du plan orbit
	@variable(modelf, -1 <= sin_meanAnomaly[1:NUM_SATELLITE] <= 1)
	@variable(modelf, -1 <= cos_meanAnomaly[1:NUM_SATELLITE] <= 1)

	# Identité trigonométrique
	@constraint(modelf, [i=1:NUM_SATELLITE], sin_inclinaison[i]^2 + cos_inclinaison[i]^2 == 1)
	@constraint(modelf, [i=1:NUM_SATELLITE], sin_noeudAscendant[i]^2 + cos_noeudAscendant[i]^2 == 1)
	@constraint(modelf, [i=1:NUM_SATELLITE], sin_meanAnomaly[i]^2 + cos_meanAnomaly[i]^2 == 1)

	# Compute the geocentric distance.
	@constraint(modelf, [i=1:NUM_SATELLITE], sum(activation_altitude[i,j] for j = 1:length(altitudeSet))==1)
	@constraint(modelf, [i=1:NUM_SATELLITE], altitude[i] == RAYON + sum(altitudeSet[j]*activation_altitude[i,j] for j=1:length(altitudeSet)))

	# Thetamax permet de déterminer la surface couverte par le satellite à partir du demi angle d'ouverture des capteurs, de l'altitude et de l'excentricité
	@constraint(modelf, [i=1:NUM_SATELLITE], Thetamax[i] == (-alphaHalf + sum(activation_altitude[i,j]*asin(((RAYON+altitudeSet[j])/RAYON)*sin(alphaHalf))
	for j=1:length(altitudeSet))))

	@NLexpression(modelf, LatitudeSat[i=1:NUM_SATELLITE, p=1:NUM_TIME], ((asin(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
	cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))*sqrt(altitude[i]/μ))))))

	@NLexpression(modelf, LongitudeSat[i=1:NUM_SATELLITE, p=1:NUM_TIME], ((((((-(angle0 + (we*((p-1)*dt))) + 
	atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
	+ (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
	+ ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
	+ (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
	((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
	- (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
	cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
	+ ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
	cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) %(2*pi)) + (2*pi))%(2*pi)))))

	@NLconstraint(modelf, [p=[time],j=1:NUM_PIX,i=1:NUM_SATELLITE], abs(- (asin(((sin_inclinaison[i]*altitude[i]*sin_meanAnomaly[i])*
	cos(sqrt(μ/(altitude[i]^3))*((p-1)*dt))/altitude[i]) + ((sin_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i])*sin(sqrt(μ/(altitude[i]^3))*((p-1)*dt))*sqrt(altitude[i]/μ)))) + Latitude) <= Thetamax[i])

	@NLconstraint(modelf, [p=[time],j=1:NUM_PIX,i=1:NUM_SATELLITE], abs(- (((((-(angle0 + (we*((p-1)*dt))) + 
	atan(((sin_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
	+ (cos_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
	+ ((-(sin_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])
	+ (cos_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))),
	((cos_noeudAscendant[i]*altitude[i]*cos_meanAnomaly[i])
	- (sin_noeudAscendant[i]*cos_inclinaison[i]*altitude[i]*sin_meanAnomaly[i]))*
	cos(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))/altitude[i])
	+ ((-(cos_noeudAscendant[i]*sqrt(μ/(altitude[i]))*sin_meanAnomaly[i])-(sin_noeudAscendant[i]*cos_inclinaison[i]*sqrt(μ/(altitude[i]))*
	cos_meanAnomaly[i]))*sin(((sqrt(μ/(altitude[i]^3)))*((p-1)*dt))*sqrt(altitude[i]/μ))))) %(2*pi)) + (2*pi))%(2*pi))) + Longitude)*cos(Latitude) <= Thetamax[i])

	@objective(modelf, Min, 0)

	status = "INFEASIBLE"
	while string(status)!="LOCALLY_SOLVED" && string(status)!="ALMOST_LOCALLY_SOLVED"

		for i=1:NUM_SATELLITE
			start_variable_inc = rand(Uniform(0,1))
			set_start_value(sin_inclinaison[i], start_variable_inc)
			set_start_value(cos_inclinaison[i], sqrt(1-start_variable_inc^2))
			start_variable_noeud = rand(Uniform(-1,1))
			set_start_value(sin_noeudAscendant[i], start_variable_noeud)
			set_start_value(cos_noeudAscendant[i], sqrt(1-start_variable_noeud^2))
			start_variable_mean = rand(Uniform(-1,1))
			set_start_value(sin_meanAnomaly[i], start_variable_mean)
			set_start_value(cos_meanAnomaly[i], sqrt(1-start_variable_mean^2))
			set_start_value(altitude[i], RAYON+altitudeSet[1])
		end

		optimize!(modelf)
		status = termination_status(modelf)
		global cputime += solve_time(modelf)
	end
	println(status)
	Theta = value.(modelf[:Thetamax])
	LatitudeSat = value.(modelf[:LatitudeSat])
	LongitudeSat = value.(modelf[:LongitudeSat])
	if Theta[1] < abs(LatitudeSat[1,time]-Latitude) || Theta[1] < abs(LongitudeSat[1,time]-Longitude)*cos(Latitude)
		println("PROBLEM")
	end
	return Theta, LatitudeSat, LongitudeSat, cputime
end

#################### LIST of PARAMETERS ######################
ObjectiveThreshold = 0.000
NUM_PIXEL = 2
SizeAddOrbitals = 1#NUM_PIXEL*n_0/4
SizeAddOrbitalsCG = NUM_PIXEL*n_0
#Latitude = [80, 61, 33, 22, 16, 8, 1, -13, -15, -24, -31, -38, -53]*pi/180
#Longitude = [358, 116, 139, 343, 340, 55, 220, 24, 100, 9, 152, 144, 46]*pi/180 # !!!!!!! Longitude ALWAYS between [0,2pi] !!!!!!!!!!!!!!!
Latitude = [32, 35, -54]*pi/180
Longitude = [42, 134, 23]*pi/180
coordsCibles = [SA[Latitude[i], Longitude[i]] for i in 1:NUM_PIXEL]
inst = DonneesSCDP(coordsCibles,alphaHalf,5,SatelliteToolbox.DateTime(1970,1,1,0),1,n_0,dt)
ecrireInst(inst)
println(Theta*180/pi)
global flag=1

#################### LIST of PARAMETERS ######################

global percentage_orbits = 0.1
global absolute_orbits_CG = 100
NbRun=1#00#0
storeResult = zeros(NbRun)
storeTime = ones(NbRun)
global size_constellation = 0
for iter in 1:NbRun
	println("iter: ",iter)
	Random.seed!(iter)

	for j in 1:NUM_PIXEL
		println("NUM_PIXEL : ", j)
		for number_period in 1:n_0
			global modelf = Model(Ipopt.Optimizer)
			set_optimizer_attribute(modelf, "print_level", 0)
			#set_optimizer_attribute(modelf, "max_iter", 10000)
			global Theta, lat_Sat, long_Sat, cputime = InitialSolution(modelf,Latitude[j],Longitude[j],1+(number_period-1)*NUM_PERIOD)
			global EDGETemp = zeros(1,n_0,NUM_PIXEL)
			observability = 0
	
			for j in 1:NUM_PIXEL
				for k in 1:n_0
					for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k)
						if Theta[1] <= abs(lat_Sat[1,p]-Latitude[j]) || Theta[1] <= abs(long_Sat[1,p]-Longitude[j])*cos(Latitude[j])
							observability += 1
						end
					end
					if observability > 0
						global EDGETemp[1,k,j] = 1
					end
				end
			end
	
			if j == 1 && number_period == 1
				global EDGE = EDGETemp
			else
				global EDGE = [EDGE;EDGETemp]
			end
	
			global CoverageSatLatTemp = zeros(1,NUM_TIME,NUM_PIXEL)
			global CoverageSatLongTemp = zeros(1,NUM_TIME,NUM_PIXEL)
			global ThetamaxTemp = zeros(1,NUM_TIME,NUM_PIXEL)
	
			for p in 1:NUM_TIME
				for t in 1:NUM_PIXEL
					global CoverageSatLatTemp[1,p,t] = abs(Latitude[t] - lat_Sat[1,p])
					global CoverageSatLongTemp[1,p,t] = abs(Longitude[t] - long_Sat[1,p])*cos(Latitude[t])
					global ThetamaxTemp[1,p,t] = Theta[1]
	
					if p == 1+(number_period-1)*NUM_PERIOD && t == j
						if CoverageSatLatTemp[1,p,t] > ThetamaxTemp[1,p,t] || CoverageSatLongTemp[1,p,t] > ThetamaxTemp[1,p,t]
							println("ERROR")
						end
					end
				end
			end
	
			if j == 1 && number_period == 1
				global CoverageSatLat = CoverageSatLatTemp
				global CoverageSatLong = CoverageSatLongTemp
				global Thetamax = ThetamaxTemp
			else
				global CoverageSatLat = [CoverageSatLat;CoverageSatLatTemp]
				global CoverageSatLong = [CoverageSatLong;CoverageSatLongTemp]
				global Thetamax = [Thetamax;ThetamaxTemp]
			end
		end
	end	

	global INDEX = (NbrTOTALsat+1):(NbrTOTALsat+n_0*NUM_PIXEL)
	global INDEX_PRICING = []
	#global CoverageSatLat,CoverageSatLong,Thetamax = getCoverageSat(INDEX,NUM_PIXEL,Latitude[1:NUM_PIXEL],Longitude[1:NUM_PIXEL])

	###################### START COLUMN GENERATION ##############
	global resultTime = @elapsed while flag>0
		#=
		println("PERCENTAGE ORBITS: ",size(INDEX_PRICING)[1]*100/NbrTOTALsat)

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

		if flag_observed == 1.0
		=#
		INDEX = collect(INDEX)
			#if size(INDEX_PRICING)[1]*100/NbrTOTALsat >= percentage_orbits
			if size(INDEX_PRICING)[1] >= absolute_orbits_CG
				println("NUM_ACTIVES_SATS : ", size(INDEX)[1])
				global model1 = Model(CPLEX.Optimizer)
				#global CoverageSatLat,CoverageSatLong,Thetamax = getCoverageSat(INDEX,NUM_PIXEL,Latitude[1:NUM_PIXEL],Longitude[1:NUM_PIXEL])
				global model1,INDEX_Original,NUM_OBS,XiTempModel1 = Original(model1,INDEX,Thetamax,CoverageSatLat,CoverageSatLong)
				INDEX = collect(INDEX)
				for i in setdiff(INDEX,INDEX_Original)
					#println(i)
					k = indexin(i,INDEX)
					deleteat!(INDEX,findall(x->x==i,INDEX))
					global EDGE = EDGE[1:end .!= k,:,:]
					global Thetamax = Thetamax[1:end .!= k,:,:]
					global CoverageSatLat = CoverageSatLat[1:end .!= k,:,:]
					global CoverageSatLong = CoverageSatLong[1:end .!= k,:,:]
				end
				println("size final constellation : ", size(INDEX)[1])
				#println(size(Thetamax))
				global storeTime[iter] += solve_time(model1)
			
				global INDEX_FINAL = INDEX
				storeResult[iter] = size(INDEX_FINAL)[1]
				break
			end
			
			if size(INDEX)[1] == 1
				global INDEX_FINAL = INDEX
				storeResult[iter] = size(INDEX_FINAL)[1]
				break
			end
			#=
			if size(INDEX)[1] == size_constellation
				global INDEX_FINAL = INDEX
				storeResult[iter] = size(INDEX_FINAL)[1]
				break
			end
			global size_constellation = size(INDEX)[1]
			=#

#			if size(INDEX)[1] == 2.0
#				global INDEX_FINAL = INDEX
#				storeResult[iter] = size(INDEX_PRICING)[1]*100/NbrTOTALsat
#				break
#			end

			# Solving Restricted Master Problem
			global model = Model(CPLEX.Optimizer)
			global model,u = RMS(model,INDEX,EDGE)
			global storeTime[iter] += solve_time(model)
			println(u)

			# Solution Pricing Problem
			global flag = 0
			addedOrbitals = 0
			for i in shuffle(1:NbrTOTALsat)
				if i ∉ INDEX_PRICING
					global CoverageSatLatTemp,CoverageSatLongTemp,ThetamaxTemp = getCoverageSat(i,NUM_PIXEL,Latitude[1:NUM_PIXEL],Longitude[1:NUM_PIXEL])
					global modelp = Model(CPLEX.Optimizer)
					set_optimizer_attribute(modelp, "CPX_PARAM_SCRIND", 0)
					# Convert INDEX vector in index of Keplerian's parameters
					IndexAltitude,Indexinclination,IndexRaan,IndexMeanAnom = ExtractKepParamFromINDEX(i)
					# Show a binary array for each target and time step to see if solution is valide
					TimeSeen_track = Validation(altitudeSet,InclinaisonSet,RAANSet,MeanAnomalySet,IndexMeanAnom,IndexRaan,Indexinclination,IndexAltitude,inst)
					global EDGETemp = zeros(1,n_0,NUM_PIXEL)
					observability = 0
					for j in 1:NUM_PIXEL
						for k in 1:n_0
							for p in (k-1)*NUM_PERIOD+1:(NUM_PERIOD*k)
								observability += TimeSeen_track[j,p]
							end
							if observability > 0
								global EDGETemp[1,k,j] = 1
							end
						end
					end
					global modelp = Pricing(modelp,ThetamaxTemp,CoverageSatLatTemp,CoverageSatLongTemp,u,EDGETemp)
					global storeTime[iter] += solve_time(modelp)
					#println(i)
					#println(termination_status(modelp))
					if string(termination_status(modelp))=="OPTIMAL"
						if objective_value(modelp) < ObjectiveThreshold
							println("objective: ", objective_value(modelp))
							push!(INDEX, i)
							push!(INDEX_PRICING, i)
							global EDGE = [EDGE;EDGETemp]
							global CoverageSatLat = [CoverageSatLat;CoverageSatLatTemp]
							global CoverageSatLong = [CoverageSatLong;CoverageSatLongTemp]
							global Thetamax = [Thetamax;ThetamaxTemp]
							global flag = 1
							addedOrbitals += 1
							if addedOrbitals >= SizeAddOrbitalsCG
								break
							end
						end
					end
				end
			end
		end
	end
#end

# Print Results
if !isdir("Results/Global/"*string(T_terre/(n_0*3600))*"H_RevisitTime/")
    mkdir("Results/Global/"*string(T_terre/(n_0*3600))*"H_RevisitTime/")
end
global name_file = "Results/Global/"*string((T_terre/(n_0*3600)))*"H_RevisitTime/ObjectiveThreshold"*string(ObjectiveThreshold)*
                "_SizeAddOrbitals"*string(SizeAddOrbitalsCG)*"_Nb_Cibles"*string(NUM_PIXEL)*"_TimeStep"*string(dt)*"_percentage"*string(percentage_orbits)*"-v3-iter.txt"
println(name_file)
open(name_file, "w") do Output
    write(Output," lat : ",  string(inst.coordsCibles[1][1]*(180/pi)))
    write(Output," ;long : ",  string(inst.coordsCibles[1][2]*(180/pi)))
    write(Output,"; ;alpha : ",string(alphaHalf*180/pi))
    write(Output,"\n"," ")
    minValue = minimum(storeResult)
    c = counter(storeResult)
    write(Output,"\n Min number of satellites : ", string(minValue))
    write(Output,";\n Rate minValue (%) : ", string(c[minValue]*100/NbRun))
	write(Output,";\n MeanValue : ", string(sum(storeResult)/NbRun))
	write(Output,";\n MeanTime : ", string(sum(storeTime)/NbRun))
    write(Output,";\n StoreResult : ", string(storeResult))
	write(Output,";\n StoreTime : ", string(storeTime))
end
