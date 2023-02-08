using EquivalentCircuits, PlotlyJS, Printf, RobustModels
using Optimization, OptimizationCMAEvolutionStrategy, OptimizationEvolutionary
# --- project pages of Evolutionary.jl, CMAEvolutionStrategy.jl: -------------------------------------------------------------------------------
# --- https://docs.sciml.ai/Optimization/stable/optimization_packages/evolutionary/
# --- https://docs.sciml.ai/Optimization/stable/optimization_packages/cmaevolutionstrategy/
# --- ----------------------------------------------------------------------------------------------------------------------
# --- Selection of optimisation methods from BlackBoxOptim package (see below): --------------------------------------------
# --- remark: no. 14 may crash "simultaneous_perturbation_stochastic_approximation"
idx_methode_selection = [4]
# --- select part of data points for equivalent circuit fitting:
idx_data = range(8, stop= 52)

# --- measurement frequencies: --------------------------------------------------------------------------------------------
frequ_data = [0.0199553, 0.0251206, 0.0316296, 0.0398258, 0.0501337, 0.0631739, 0.0794492, 0.1001603, 0.1260081, 0.1588983, 
    0.2003205, 0.2520161, 0.316723, 0.400641, 0.5040323, 0.6334459, 0.7923428, 0.999041, 1.266892, 1.584686, 1.998082, 
    2.504006, 3.158693, 3.945707, 5.008013, 6.317385, 7.944915, 9.93114, 12.40079, 15.625, 19.86229, 24.93351, 31.25, 
    38.42213, 50.22321, 63.3446, 79.00281, 100.4464, 125.558, 158.3615, 198.6229, 252.4038, 315.5048, 397.9953, 505.5147, 
    627.7902, 796.875, 998.264, 1265.625, 1577.524, 1976.103, 2527.573]
# --- measured values impedance values (normed by max value): ------------------------------------------------------------
Z_data = ComplexF64[1.0 + 0.00033364003270441563im, 0.9959327788858627 - 0.004310226899902616im, 0.9942611041385668 - 0.005671195638022529im, 
    0.9896626966561299 - 0.00856668211620483im, 0.9848663962129538 - 0.011108044348854046im, 0.9839133853757104 - 0.011680892393098745im, 
    0.9803825583393658 - 0.014357655072569432im, 0.9780390890674556 - 0.016721955182452106im, 0.977945350296579 - 0.019502872051785466im, 
    0.9749352942095479 - 0.023195138082417213im, 0.973617743707785 - 0.027590444894622im, 0.9703941715315353 - 0.03213677528212766im, 
    0.9655562094124143 - 0.03752154689802784im, 0.9597964827130084 - 0.04379162912776074im, 0.9520005416017874 - 0.04925451638605793im, 
    0.9442670930044841 - 0.054795519286752116im, 0.9348307234029257 - 0.05974284330522906im, 0.9245767434109456 - 0.06472662129015795im, 
    0.9139634314639392 - 0.06913234352134902im, 0.9049228477838593 - 0.07387656686959375im, 0.8949136301380565 - 0.07923529993802828im, 
    0.8851127208719789 - 0.08469818719632546im, 0.8739630148471799 - 0.0921400040620134im, 0.8622300453591497 - 0.1000505147820834im, 
    0.8483046301745104 - 0.10995557823802357im, 0.8333168422532719 - 0.1212302692906579im, 0.8155116834962478 - 0.13342151721408374im, 
    0.7950141389312739 - 0.14564921910396153im, 0.7706472662129017 - 0.15724678814517012im, 0.7418694635538452 - 0.16734453685235623im, 
    0.708347437546544 - 0.1744582680199768im, 0.6755492831587884 - 0.17673924477796932im, 0.6430635913406207 - 0.17507798544965972im, 
    0.6145565895752072 - 0.17028168500648364im, 0.580763762674263 - 0.16053285283533744im, 0.5548293693984576 - 0.1501643032345084im, 
    0.5330351051696932 - 0.13955099128750204im, 0.5123240445155008 - 0.12800029163173163im, 0.4953156653109263 - 0.11759528806445062im, 
    0.4794894361612932 - 0.10735172349145677im, 0.4655952672336127 - 0.09793097701837802im, 0.45218020757930044 - 0.08854147680225807im, 
    0.44080136233680345 - 0.08018310306577858im, 0.4300318191050031 - 0.07178306765335403im, 0.4198247084984612 - 0.06336220140295694im, 
    0.41134134973414643 - 0.055706868448050506im, 0.40291527577425623 - 0.04717664029829761im, 0.39567135186930735 - 0.03873494320992798im, 
    0.3888857063997542 - 0.029465220312150108im, 0.38332908037058067 - 0.020200705123865372im, 0.37827760216224104 - 0.009806116975570635im, 
    0.3737000255177766 + 0.002913495258380507im]

# --- reference fit is computed via a commercial EIS-software: --------------------------------------------------------------
trace_ref_fit_SW = ComplexF64[0.9788651042462025 - 0.003427725046892952im, 0.9787042129564272 - 0.0041986283934134444im, 0.9785003527527543 - 0.005145223347737662im, 
    0.9782393462484965 - 0.006306332111776151im, 0.9779013997203362 - 0.007727852770328056im, 0.9774552153205819 - 0.009476179242838563im, 
    0.9768659237035731 - 0.011593226994444013im, 0.9760613725326219 - 0.01419779962088405im, 0.9749709741261627 - 0.017323133776195227im, 
    0.9734501264062743 - 0.02111374980931685im, 0.9713406718464188 - 0.025603088443209795im, 0.9684505314612848 - 0.030769668472705174im, 
    0.9645292196704338 - 0.036566648193514566im, 0.9591448671132816 - 0.043022435797319084im, 0.9523887936368954 - 0.04948223664885959im, 
    0.9442179647804041 - 0.05564232645454433im, 0.9350958290630905 - 0.06105262795583616im, 0.924987982855181 - 0.06586070781401777im, 
    0.9145015589339538 - 0.0701685599788112im, 0.9048243526644042 - 0.07415006287794139im, 0.8950106684802762 - 0.07881112049424521im, 
    0.8853693610889181 - 0.08439055008090053im, 0.8748621766817914 - 0.09155553424484979im, 0.8636873268852573 - 0.09990018322926661im, 
    0.8498691726172095 - 0.11037443576081168im, 0.8339703638271251 - 0.12180969637711742im, 0.8154128823748253 - 0.13378783963504204im, 
    0.7942784123186765 - 0.1454474185257693im, 0.7701271385576475 - 0.1562531797843441im, 0.7418846967259558 - 0.16570911888871886im, 
    0.7098260570782756 - 0.17258871607415407im, 0.6778320548678061 - 0.1755989237878614im, 0.6456734302156821 - 0.17484887477148447im, 
    0.6169363253992103 - 0.17095547213511444im, 0.5820365788624768 - 0.16196409353188074im, 0.5548571975533133 - 0.1514921236560258im, 
    0.5320564866437297 - 0.1401634302966748im, 0.5107034591263453 - 0.12731820019180443im, 0.49382330842857997 - 0.11564715504857762im, 
    0.4788316298662813 - 0.10435619961101911im, 0.4661622167087829 - 0.09449953735021457im, 0.4541975586774136 - 0.08541745893062654im, 
    0.44379749463126544 - 0.07803589286507319im, 0.43325665751158693 - 0.07102760719786393im, 0.4224787798478673 - 0.06387538861048644im, 
    0.41289733172108517 - 0.05686348014163662im, 0.4029795762886971 - 0.04810976004278965im, 0.3947009737016643 - 0.03873887736795084im, 
    0.3874833173425731 - 0.027874091160924586im, 0.38222507364956054 - 0.01707695883305848im, 0.3781417585287105 - 0.005417302648024877im, 
    0.37489234021412426 + 0.008160231785984899im]

# --- defined equavalent circuit model (notation according to Julia Package: "EquivalentCircuits")
circ_strg_ref = "R1-L2-[P3,R4]-[P5,R6]-[P7,R8]"

# --- Sample of available optimisation methods from BlackBoxOptim package:
ES_methods_all =  ["01#Genetic_Algorithm_optimizer(GA)", "02#Differential_Evolution_optimizer(DE)", "03#Evolution_Strategy_algorithm(ES)", 
"04#Covariance_Matrix_Adaptation(CMAES)"]

ES_methods_sel = ES_methods_all[idx_methode_selection]

# --- function "denumber_circuit()" ----------------------------------------------------------------------------------------
denumber_circuit(Circuit) = replace(Circuit, r"[0-9]" => "")

# --- Function used to optimise the circuit parameters, using a Evolution Strategy via package "Optimization.jl": ----------
function parameopt(_circuitstring::String, _Z_measured, _frequencies, _optim_method)
    println("DBG parameopt, Begin: ", '-'^100)
    elements            = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""], init = denumber_circuit(_circuitstring))
    initial_parameters  = EquivalentCircuits.flatten(EquivalentCircuits.karva_parameters(elements));
    circfunc            = circuitfunction(_circuitstring)
    objective           = EquivalentCircuits.objectivefunction(circfunc, _Z_measured, _frequencies) 
    lower               = zeros(length(initial_parameters))
    upper               = get_parameter_upper_bound(_circuitstring)

    ### First step ###
    SR = Array{Tuple{Float64, Float64}, 1}(undef, length(initial_parameters))
    for (e, (l, u)) in enumerate(zip(lower, upper))
        SR[e] = (l, u)
    end
    # --- DBG:
    println("DBG parameopt, SR", SR)
    println("DBG parameopt, lower", lower)
    println("DBG parameopt, upper", upper)
    println("DBG parameopt, typeof(objective): ", typeof(objective))
    # ---
    if cmp(_optim_method, "04#Covariance_Matrix_Adaptation(CMAES)") == 0
        rosenbrock(x, p) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
        x0 = zeros(2)
        p = [1.0, 100.0]
        f = OptimizationFunction(rosenbrock)
        prob = Optimization.OptimizationProblem(f, x0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
        sol = solve(prob, Evolutionary.CMAES(μ = 40, λ = 100))
        # ???
        # opt_func = OptimizationFunction(objective)
        # println("DBG parameopt, typeof(opt_func)", typeof(opt_func))

        # prob = Optimization.OptimizationProblem(objective, initial_parameters, SciMLBase.NullParameters(),; lb = lower, ub = upper)
        # sol = solve(prob, Evolutionary.CMAES(μ = 40, λ = 100))
        println("DBG parameopt, sol.u: ", sol.u, "\n --- END ---", '-'^110, "\n") 
        return parameters = initial_parameters
    else
        error("func parameopt: Methode not yet implemented!")
    end
    
    # res = bboptimize(objective; SearchRange = SR, Method = _optim_method, TraceMode = :silent); #MaxSteps=70000,
    # initial_parameters  = best_candidate(res);
    # if _optim_method == :simultaneous_perturbation_stochastic_approximation
    #     initial_parameters = clamp.(initial_parameters, 0.0, 1.0)
    # end
    # fitness_1           = best_fitness(res);
    # ### Second step ###
    # inner_optimizer     = NelderMead()
    # println("DBG: [min..max] of initial_parameters: [", minimum(initial_parameters), " .. ", maximum(initial_parameters), "], _optim_method: :", _optim_method)
    # results             = optimize(objective, lower, upper, initial_parameters, Fminbox(inner_optimizer), 
    #                         Optim.Options(time_limit = 50.0)); #20.0
    # fitness_2           = results.minimum
    # best                = results.minimizer
    # parameters          = fitness_2 < fitness_1 ? best : initial_parameters
    # ---
    return parameters
end

# Optional: The parameter bounds for this particular application can be made more specific than the bounds 
# that are generally used in the package, which is application-agnostic.
function get_parameter_upper_bound(readablecircuit::String)
    elements = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""], init = denumber_circuit(readablecircuit))
    ranges = Dict('R'=>400,'C'=>0.01,'L'=>5,'P'=>[400, 1],'W'=>400,'+'=>0,'-'=>0) 
    return EquivalentCircuits.flatten([ranges[e] for e in elements])
end

# Optimise the circuit parameters and simulate the resulting impedance spectrum.
function optimise_and_simulate(_circuitstring, _Z_measured, _frequencies, _optim_method)
    _parameters  = parameopt(_circuitstring, _Z_measured, _frequencies, _optim_method)
    _trace       = simulateimpedance_noiseless(circuitfunction(_circuitstring), _parameters, _frequencies)
    return _trace, _parameters
end

#Simulate impedance spectra for circuit parameters using the considered optimisation methods.
function generate_traces(_circuitstring, _ref_data, _Z_measured, _frequencies, _opt_methods)
    _names      = Vector{String}(undef, length(_opt_methods) + 2)
    _fited_par  = Dict{Int, Union{Missing, Any}}()
    _traces     = Vector{Vector{ComplexF64}}(undef, length(_opt_methods) + 2)
    _traces[1]  = _Z_measured;   _names[1]   = "Z_measured";       _fited_par[1] = missing  
    _traces[2]  = _ref_data;     _names[2]   = "Ref_SW_fit";       _fited_par[2] = missing
    for _i in 3:length(_traces)
        @info(string("--- generate_traces, #: ", _i - 2, ", next methode: ", _opt_methods[_i - 2], "\t-----------------------------"))
        _traces[_i], _fited_par[_i] = optimise_and_simulate(_circuitstring, _Z_measured, _frequencies, _opt_methods[_i - 2])
        _names[_i]  = _opt_methods[_i - 2]
    end
    return _traces, _names, _fited_par
end

# The fitting errors calculated using the modulus weighted objective function,
# you can adjust the function to see other fitting quality metrics (e.g. removal of the denominator gives the MSE).
function trace_quality(_Z_measured::Vector{ComplexF64}, trace::Vector{ComplexF64}; 
    weighing_fac::Vector{Float64}=Vector{Float64}(undef, 0), b_tail::Bool=true)
    _n_measuremnt = length(_Z_measured); _n_weighing = length(weighing_fac)
    if isempty(weighing_fac) 
        weighing_fac = ones(_n_measuremnt, )
    elseif _n_weighing  < _n_measuremnt
        fill_vec = ones(_n_measuremnt - _n_weighing, )
        if b_tail
            weighing_fac = vcat(fill_vec, weighing_fac)
        else
            weighing_fac = vcat(weighing_fac, fill_vec)
        end
    elseif length(weighing_fac)  > _n_measuremnt
        @warn("Weighing factor is larger then number of measured points! - Truncated weighing vector will be used!")
        weighing_fac = weighing_fac[1:_n_measuremnt]
    end
    return mean((abs.(weighing_fac .* (_Z_measured - trace)).^2)./(abs.(_Z_measured).^2 .+ abs.(trace).^2))
end

# --- Nyquist plots of arbitraty number of impedance spectra, compared to Z_measured and reference fit.
function plot_Nyquist(_traces::Vector{Vector{ComplexF64}}, names::Vector{String}, _frequ_data::Vector{Float64}) 
    _dtick = 1.0e-1
    s_pts_info = []
    for i_ndx in eachindex(_frequ_data)
        push!(s_pts_info, @sprintf("#:%i, f=%.2fHz", i_ndx, _frequ_data[i_ndx]))
    end
    # --- traces:
    T = Array{GenericTrace{Dict{Symbol, Any}}}(undef,length(_traces))
    for i in 1:length(_traces)
        if i == 1
        T[i] = PlotlyJS.scatter(; x= real(_traces[i]),    y=  imag(_traces[i]),    name = names[i], text = s_pts_info, mode = "markers")
        else
        T[i] = PlotlyJS.scatter(; x= real(_traces[i]),    y=  imag(_traces[i]),    name = names[i])
        end
    end
    # --- layout:
    plt_layout = PlotlyJS.Layout(;
        title_text          = "Comparison of Fit Results",
        xaxis_title_text    = "normed z<sub>Real</sub> / -",
        xaxis_dtick         = _dtick,
        yaxis_title_text    = "normed z<sub>Imag</sub> / -",
        yaxis_dtick         = _dtick,
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
    )
    return PlotlyJS.Plot(T, plt_layout)
end

# --- Plot the fitting errors to compare the different optimisation methods. -----------------------------------------------
function Plot_fitting_errors(_method_names, fit_errors)
    bar_plot = PlotlyJS.bar(x= _method_names, y= fit_errors, text= string.(fit_errors), textposition= "auto")
    plt_layout = PlotlyJS.Layout(;
        title_text          = "Comparison of fitting methods",
        yaxis_title_text    = "Fitting error",
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
    )
    return PlotlyJS.Plot(bar_plot,plt_layout)
end
# --- main part: -----------------------------------------------------------------------------------------------------------
n_measured = length(frequ_data)
if idx_data[end] > n_measured
    error("Selection of measured points \"idx_data\" exceeds vector of measured points!")
end
frequ_data_sel       = frequ_data[idx_data]
Z_data_sel           = Z_data[idx_data]
trace_ref_fit_SW_sel = trace_ref_fit_SW[idx_data]
# --- Evaluate the optimisation methods and make plots: --------------------------------------------------------------------
@time begin
    traces_, names_, fitted_params = generate_traces(circ_strg_ref, trace_ref_fit_SW_sel, Z_data_sel, frequ_data_sel, ES_methods_sel);
end
println("--- END generate traces: ", '-'^100)
# --- list qualities / deviation from measured values and sort fit results according to their quality: ---------------------
# Fitting_err         = [trace_quality(Z_data_sel, traces_[i]) for i in 1:length(traces_)]
# idx_sorted          = sortperm(Fitting_err[3:end]); idx_sorted = 2 .+ idx_sorted
# Fitting_err_sorted  = vcat(Fitting_err[1:2], Fitting_err[idx_sorted])
# names_sorted        = vcat(names_[1:2], names_[idx_sorted])
# traces_sorted       = vcat(traces_[1:2], traces_[idx_sorted])
# # --- plot:
# quality_plot        = Plot_fitting_errors(names_sorted, Fitting_err_sorted)
# fn_plt = raw"c:\tmp\plt\quality_plot.html"
# PlotlyJS.savefig(quality_plot, fn_plt)

# nyqs_plot           = plot_Nyquist(traces_sorted[1:4], names_sorted[1:4], frequ_data_sel)
# fn_plt = raw"c:\tmp\plt\nyquist_plot.html"
# PlotlyJS.savefig(nyqs_plot, fn_plt)

