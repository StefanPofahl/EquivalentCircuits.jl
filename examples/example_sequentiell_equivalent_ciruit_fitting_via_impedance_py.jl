using PyCall, EquivalentCircuits, RobustModels, Printf, PlotlyJS
using StefanPHyCentaImpedance # MyLibOptimizeEquivalentCircuit

# --- measured values:
frequ_data = [0.0199553, 0.0251206, 0.0316296, 0.0398258, 0.0501337, 0.0631739, 0.0794492, 0.1001603, 0.1260081, 0.1588983, 
    0.2003205, 0.2520161, 0.316723, 0.400641, 0.5040323, 0.6334459, 0.7923428, 0.999041, 1.266892, 1.584686, 1.998082, 
    2.504006, 3.158693, 3.945707, 5.008013, 6.317385, 7.944915, 9.93114, 12.40079, 15.625, 19.86229, 24.93351, 31.25, 
    38.42213, 50.22321, 63.3446, 79.00281, 100.4464, 125.558, 158.3615, 198.6229, 252.4038, 315.5048, 397.9953, 505.5147, 
    627.7902, 796.875, 998.264, 1265.625, 1577.524, 1976.103, 2527.573]
Z_data = ComplexF64[0.0192023 + 6.406656e-6im, 0.0191242 - 8.276627e-5im, 0.0190921 - 0.0001089im, 0.0190038 - 0.0001645im, 
    0.0189117 - 0.0002133im, 0.0188934 - 0.0002243im, 0.0188256 - 0.0002757im, 0.0187806 - 0.0003211im, 0.0187788 - 0.0003745im, 
    0.018721 - 0.0004454im, 0.0186957 - 0.0005298im, 0.0186338 - 0.0006171im, 0.0185409 - 0.0007205im, 0.0184303 - 0.0008409im, 
    0.0182806 - 0.0009458im, 0.0181321 - 0.0010522im, 0.0179509 - 0.0011472im, 0.017754 - 0.0012429im, 0.0175502 - 0.0013275im, 
    0.0173766 - 0.0014186im, 0.0171844 - 0.0015215im, 0.0169962 - 0.0016264im, 0.0167821 - 0.0017693im, 0.0165568 - 0.0019212im, 
    0.0162894 - 0.0021114im, 0.0160016 - 0.0023279im, 0.0156597 - 0.002562im, 0.0152661 - 0.0027968im, 0.0147982 - 0.0030195im, 
    0.0142456 - 0.0032134im, 0.0136019 - 0.00335im, 0.0129721 - 0.0033938im, 0.0123483 - 0.0033619im, 0.0118009 - 0.0032698im, 
    0.011152 - 0.0030826im, 0.010654 - 0.0028835im, 0.0102355 - 0.0026797im, 0.0098378 - 0.0024579im, 0.0095112 - 0.0022581im, 
    0.0092073 - 0.0020614im, 0.0089405 - 0.0018805im, 0.0086829 - 0.0017002im, 0.0084644 - 0.0015397im, 0.0082576 - 0.0013784im, 
    0.0080616 - 0.0012167im, 0.0078987 - 0.0010697im, 0.0077369 - 0.0009059im, 0.0075978 - 0.0007438im, 0.0074675 - 0.0005658im, 
    0.0073608 - 0.0003879im, 0.0072638 - 0.0001883im, 0.0071759 + 5.594581e-5im]  

# --- generate frequency vector with n_elements with the same range_data as given in the measurement:
n_elements = 1000
# frequ_vec_all   = exp10.(LinRange(log10(frequ_data_all[1]), log10(frequ_data_all[end]), n_elements))
frequ_simu       = exp10.(LinRange(log10(frequ_data[1]), log10(frequ_data[end]), n_elements))

# --- Optimized Parameters Calculated by GAMRY-Software: -------------------------------------------------------------------
# --- Reference case: ------------------------------------------------------------------------------------------------------
circ_strg_ref = "R1-L2-[P3,R4]-[P5,R6]-[P7,R8]"
const R1_ref    = 0.007031                                          # Gamry: HFR
const L2_ref    = 0.00000004257                                     # Gamry: Lstray
# ---
const P3w_ref   = 149.9;            const P3n_ref   = 0.9763        # Gamry: Q_3, a_3
const R4_ref    = 0.00132                                           # Gamry: R_3
# ---
const P5w_ref   = 1.948;            const P5n_ref   = 0.7817        # Gamry: Q_2, a_2
const R6_ref    = 0.009341                                          # Gamry: R2
# ---
const P7w_ref   = 0.2224;           const P7n_ref   = 9.97E-01      # Gamry: Q_1, a_1
const R8_ref    = 0.001118                                          # Gamry: R_1

params_Gamry = (R1 = R1_ref, L2 = L2_ref, P3w = P3w_ref, P3n = P3n_ref, R4 = R4_ref, 
    P5w = P5w_ref, P5n = P5n_ref, R6 = R6_ref, P7w = P7w_ref, P7n = P7n_ref, R8 = R8_ref)    
# --- END Reference case: --------------------------------------------------------------------------------------------------

# --- local functions: -----------------------------------------------------------------------------------------------------
function fit_three_RP_EC(_circ_strg_in::AbstractString, _eq_params_in::NamedTuple, _frequ::Vector{Float64}, _Z_data::Vector{ComplexF64}, b_DBG::Bool=false)
    _circuits = PyCall.pyimport("impedance.models.circuits")
    _circuit_ImpPy  = MylibExpECircJLStrToImpPy(_circ_strg_in)
    _Q_opt = +Inf; _Q = +Inf
    fact_ = 1.3 # facto to modify initial guess parameters
    n_optimization_loops = 1 # number of for-loops to optimize parameters
    # --- constants:
    idxR1       = 1
    idxL2       = 2
    idxCPE3_0   = 3
    idxCPE3_1   = 4
    idxR4       = 5
    idxCPE5_0   = 6
    idxCPE5_1   = 7
    idxR6       = 8
    idxCPE7_0   = 9
    idxCPE7_1   = 10
    idxR8       = 11
    # -------------------------------------------------
    function _sign()   
        if rand(Int8, 1)[] > 0
            sign_ = 1
        else
            sign_ = -1
        end
        return sign_
    end
        
#---  R1=1 L2=2 P3_0=3 P3_1=4 R4=5 P5_0=6 P5_1=7 R6=8 P7_0=9 P7_1=10 R8=11    
    function _const_all_but_R1_L2(_b::Bool=false)
        _par_opt = fact_ .* _params_ImPy[[1,2]]
        _dict_const = Dict(
        "CPE3_0" => _params_ImPy[idxCPE3_0],    "CPE3_1" => _params_ImPy[idxCPE3_1],    "R4" => _params_ImPy[idxR4], 
        "CPE5_0" => _params_ImPy[idxCPE5_0],    "CPE5_1" => _params_ImPy[idxCPE5_1],    "R6" => _params_ImPy[idxR6], 
        "CPE7_0" => _params_ImPy[idxCPE7_0],    "CPE7_1" => _params_ImPy[idxCPE7_1],    "R8" => _params_ImPy[idxR8], 
        )
        _customCircuit   = _circuits.CustomCircuit(initial_guess= _par_opt, constants= _dict_const, circuit= _circuit_ImpPy)
        if _b; println("_const_all_but_R1_L2: ", '-'^100); end
        if _b; py"print"(_customCircuit); end
        _customCircuit.fit(_frequ, _Z_data)
        if _b; py"print"(_customCircuit); end
        _Z_simulated = _customCircuit.predict(_frequ)
        return vcat(_customCircuit.parameters_, _params_ImPy[3:end]), _Z_simulated
    end
#---  R1=1 L2=2 P3_0=3 P3_1=4 R4=5 P7_0=9 P7_1=10 R8=11    
    function _const_all_but_CPE3_R4(_b::Bool=false)
        _par_opt = fact_ .* _params_ImPy[[3,4,5]]
        _dict_const =  Dict(
            "R1"     => _params_ImPy[idxR1],        "L2"     => _params_ImPy[idxL2], 
            "CPE5_0" => _params_ImPy[idxCPE5_0],    "CPE5_1" => _params_ImPy[idxCPE5_1],    "R6" => _params_ImPy[idxR6], 
            "CPE7_0" => _params_ImPy[idxCPE7_0],    "CPE7_1" => _params_ImPy[idxCPE7_1],    "R8" => _params_ImPy[idxR8], 
            )
        _customCircuit   = _circuits.CustomCircuit(initial_guess= _par_opt, constants= _dict_const, circuit= _circuit_ImpPy)
        if _b; println("_const_all_but_CPE3_R4: ", '-'^100); end
        if _b; py"print"(_customCircuit); end
        _customCircuit.fit(_frequ, _Z_data)
        if _b; py"print"(_customCircuit); end
        _Z_simulated = _customCircuit.predict(_frequ)
        return vcat(_params_ImPy[1:2], _customCircuit.parameters_, _params_ImPy[6:end]), _Z_simulated
    end
#---  R1=1 L2=2 P3_0=3 P3_1=4 R4=5 P5_0=6 P5_1=7 R6=8 P7_0=9 P7_1=10 R8=11    
    function _const_all_but_CPE5_R6(_b::Bool=false)
        _par_opt = fact_ .* _params_ImPy[[6,7,8]]
        _dict_const =   Dict(
        "R1"     => _params_ImPy[idxR1],        "L2"     => _params_ImPy[idxL2], 
        "CPE3_0" => _params_ImPy[idxCPE3_0],    "CPE3_1" => _params_ImPy[idxCPE3_1],    "R4" => _params_ImPy[idxR4], 
        "CPE7_0" => _params_ImPy[idxCPE7_0],    "CPE7_1" => _params_ImPy[idxCPE7_1],    "R8" => _params_ImPy[idxR8], 
        )
        _customCircuit   = _circuits.CustomCircuit(initial_guess= _par_opt, constants= _dict_const, circuit= _circuit_ImpPy)
        if _b; println("_const_all_but_CPE5_R6: ", '-'^100); end
        if _b; py"print"(_customCircuit); end
        _customCircuit.fit(_frequ, _Z_data)
        if _b; py"print"(_customCircuit); end
        _Z_simulated = _customCircuit.predict(_frequ)
        return vcat(_params_ImPy[1:5], _customCircuit.parameters_, _params_ImPy[9:end]), _Z_simulated
    end
#---  R1=1 L2=2 P3_0=3 P3_1=4 R4=5 P5_0=6 P5_1=7 R6=8 P7_0=9 P7_1=10 R8=11    
    function _const_all_but_CPE7_R8(_b::Bool=false)
        _par_opt = fact_ .* _params_ImPy[[9,10,11]]
        _dict_const =    Dict(
        "R1"     => _params_ImPy[idxR1],        "L2"     => _params_ImPy[idxL2], 
        "CPE3_0" => _params_ImPy[idxCPE3_0],    "CPE3_1" => _params_ImPy[idxCPE3_1],    "R4" => _params_ImPy[idxR4], 
        "CPE5_0" => _params_ImPy[idxCPE5_0],    "CPE5_1" => _params_ImPy[idxCPE5_1],    "R6" => _params_ImPy[idxR6], 
        )
        _customCircuit   = _circuits.CustomCircuit(initial_guess= _par_opt, constants= _dict_const, circuit= _circuit_ImpPy)
        if _b; println("_const_all_but_CPE7_R8: ", '-'^100); end
        if _b; py"print"(_customCircuit); end
        _customCircuit.fit(_frequ, _Z_data)
        if _b; py"print"(_customCircuit); end
        _Z_simulated = _customCircuit.predict(_frequ)
        return vcat(_params_ImPy[1:8], _customCircuit.parameters_), _Z_simulated
    end
#---  R1=1 L2=2 P3_0=3 P3_1=4 R4=5 P5_0=6 P5_1=7 R6=8 P7_0=9 P7_1=10 R8=11    
    function _const_all_but_Rx(_b::Bool=false)
        _par_opt = fact_ .* _params_ImPy[[1, 5, 8, 11]]
        _dict_const =     Dict(
        "L2" => _params_ImPy[idxL2], 
        "CPE3_0" => _params_ImPy[idxCPE3_0],    "CPE3_1" => _params_ImPy[idxCPE3_1], 
        "CPE5_0" => _params_ImPy[idxCPE5_0],    "CPE5_1" => _params_ImPy[idxCPE5_1],    
        "CPE7_0" => _params_ImPy[idxCPE7_0],    "CPE7_1" => _params_ImPy[idxCPE7_1],
        )
        _customCircuit   = _circuits.CustomCircuit(initial_guess= _par_opt, constants= _dict_const, circuit= _circuit_ImpPy)
        if _b; println("_const_all_but_Rx: ", '-'^100); end
        if _b; py"print"(_customCircuit); end
        _customCircuit.fit(_frequ, _Z_data)
        if _b; py"print"(_customCircuit); end
        _Z_simulated = _customCircuit.predict(_frequ)
        return vcat(_customCircuit.parameters_[1], _params_ImPy[2:4], _customCircuit.parameters_[2],
            _params_ImPy[6:7], _customCircuit.parameters_[3], _params_ImPy[9:end], ), _Z_simulated
    end
#---  R1=1 L2=2 P3_0=3 P3_1=4 R4=5 P5_0=6 P5_1=7 R6=8 P7_0=9 P7_1=10 R8=11    
    function _const_all_but_CPEx_0(_b::Bool=false)
        _par_opt = fact_ .* _params_ImPy[[3,6,9]]
        _dict_const =  Dict(
        "R1"     => _params_ImPy[idxR1],            "L2" => _params_ImPy[idxL2], 
        "R4"     => _params_ImPy[idxR4],            "R6" => _params_ImPy[idxR6],            "R8" => _params_ImPy[idxR8], 
        "CPE3_1" => _params_ImPy[idxCPE3_1],    "CPE5_1" => _params_ImPy[idxCPE5_1],    "CPE7_1" => _params_ImPy[idxCPE7_1], 
        )
        _customCircuit   = _circuits.CustomCircuit(initial_guess= _par_opt, constants= _dict_const, circuit= _circuit_ImpPy)
        if _b; println("_const_all_but_CPEx_0: ", '-'^100); end
        if _b; py"print"(_customCircuit); end
        _customCircuit.fit(_frequ, _Z_data)
        if _b; py"print"(_customCircuit); end
        _Z_simulated = _customCircuit.predict(_frequ)
        return vcat(_params_ImPy[1:2], _customCircuit.parameters_[1],
            _params_ImPy[4:5], _customCircuit.parameters_[2], 
            _params_ImPy[7:8], _customCircuit.parameters_[3], _params_ImPy[10:11]), _Z_simulated
    end
#---  R1=1 L2=2 P3_0=3 P3_1=4 R4=5 P5_0=6 P5_1=7 R6=8 P7_0=9 P7_1=10 R8=11    
    function _const_all_but_CPEx_1(_b::Bool=false)
        _par_opt = fact_ .* _params_ImPy[[4,7,10]]
        _dict_const =  Dict(
        "R1"     => _params_ImPy[idxR1],            "L2" => _params_ImPy[idxL2], 
        "R4"     => _params_ImPy[idxR4],            "R6" => _params_ImPy[idxR6],            "R8" => _params_ImPy[idxR8], 
        "CPE3_0" => _params_ImPy[idxCPE3_0],    "CPE5_0" => _params_ImPy[idxCPE5_0],    "CPE7_0" => _params_ImPy[idxCPE7_0], 
        )
        _customCircuit   = _circuits.CustomCircuit(initial_guess= _par_opt, constants= _dict_const, circuit= _circuit_ImpPy)
        if _b; println("_const_all_but_CPEx_1 ", '-'^100); end
        if _b; py"print"(_customCircuit); end
        _customCircuit.fit(_frequ, _Z_data)
        if _b; py"print"(_customCircuit); end
        _Z_simulated = _customCircuit.predict(_frequ)
        return vcat(_params_ImPy[1:3], _customCircuit.parameters_[1],
            _params_ImPy[5:6], _customCircuit.parameters_[2], 
            _params_ImPy[8:9], _customCircuit.parameters_[3]), _Z_simulated
    end
    
    # ---
    println("_circuit_ImpPy: ", _circuit_ImpPy) # R1-L2-p(CPE3,R4)-p(CPE5,R6)-p(CPE7,R8)
    _params_ImPy    = collect(_eq_params_in)
    _customCircuit   = _circuits.CustomCircuit(initial_guess=_params_ImPy, circuit=_circuit_ImpPy)
    if b_DBG; println("Initial fit of all parameters: ", '-'^100); end
    if b_DBG; py"print"(_customCircuit); end
    _customCircuit.fit(_frequ, _Z_data)
    if b_DBG; py"print"(_customCircuit); end
    _Z_simulated = _customCircuit.predict(_frequ)
    _Q_init = RobustModels.mean(abs.(_Z_data - _Z_simulated))
    println("_Q_init: \t\t", _Q_init)
    
    for i_ = 1:n_optimization_loops
        println("#: ", i_)
        _EQ_params, _Z_simulated = _const_all_but_R1_L2()
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_const_all_but_R1_L2: \t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end

        _EQ_params_opt, _Z_simulated = _const_all_but_CPE3_R4(true)
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_const_all_but_CPE3_R4: \t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end
        
        _EQ_params_opt, _Z_simulated = _const_all_but_CPE5_R6()
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_const_all_but_CPE5_R62: \t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end
        
        _EQ_params_opt, _Z_simulated = _const_all_but_CPE7_R8()
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_all_but_CPE7_R8:  \t\t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end
        
        _EQ_params_opt, _Z_simulated = _const_all_but_Rx()
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_all_but_Rx:       \t\t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end
        
        _EQ_params_opt, _Z_simulated = _const_all_but_CPEx_0()
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_all_but_CPEx_0:     \t\t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end
        
        _EQ_params_opt, _Z_simulated = _const_all_but_CPEx_1()
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_all_but_CPEx_1: \t\t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end

        _customCircuit   = _circuits.CustomCircuit(initial_guess= _params_ImPy, circuit= _circuit_ImpPy)
        _customCircuit.fit(_frequ, _Z_data)
        _Z_simulated = _customCircuit.predict(_frequ)
        _Q = RobustModels.mean(abs.(_Z_data - _Z_simulated))
        println("_Q_all: \t\t\t", _Q, "\t, \tΔQ: \t", _Q - _Q_opt)
        if _Q < _Q_opt
            _params_ImPy = _EQ_params
            _Q_opt = _Q
        end    

    end
    # ---
    _EQ_params_opt_NT = MyLibEqCircNamedTuple(_circ_strg_in, _params_ImPy)

    # ---
    return _Q, _EQ_params_opt_NT, _Z_simulated
end

function simulate_impedance(_circ_strg_ref::String, _ciruit_params_NT::NamedTuple,  
    _frequ_data::Vector{Float64}, _Z_data::Vector{ComplexF64})
    # ---
    _circuits = PyCall.pyimport("impedance.models.circuits")
    _circ_str_ImpPy = MylibExpECircJLStrToImpPy(_circ_strg_ref)
    _circuit_ImpPy  = _circuits.CustomCircuit(initial_guess= collect(_ciruit_params_NT), circuit= _circ_str_ImpPy)
    # --- simulate impedance according to given parameters (no optimization):
    _Z_data_simu    = _circuit_ImpPy.predict(_frequ_data, use_initial = true)
    _q_ref          = RobustModels.mean(abs.(_Z_data - _Z_data_simu))
    return _Z_data_simu, _q_ref
end

function plot_Nyquist(_Z_data::Vector{ComplexF64}, _Z_Gamry::Vector{ComplexF64}, _Z_GenAlg1::Vector{ComplexF64}, _Z_GenAlg2::Vector{ComplexF64}, _frequ_data::Vector{Float64})
    _dtick = 2.0e-3
    s_pts_info = []
    for i_ndx in eachindex(_frequ_data)
        push!(s_pts_info, @sprintf("#:%i, f=%.2fHz", i_ndx, _frequ_data[i_ndx]))
    end
    # --- traces:
    _trace_data     = PlotlyJS.scatter(; x= real(_Z_data),    y=  imag(_Z_data),    name = "measured", text = s_pts_info, mode = "markers")
    _trace_Gamry    = PlotlyJS.scatter(; x= real(_Z_Gamry),   y=  imag(_Z_Gamry),   name = "Gamry")
    _trace_Init     = PlotlyJS.scatter(; x= real(_Z_GenAlg1), y=  imag(_Z_GenAlg1), name = "Init")
    _trace_Optim    = PlotlyJS.scatter(; x= real(_Z_GenAlg2), y=  imag(_Z_GenAlg2), name = "Opt")
    # --- layout:
    plt_layout = PlotlyJS.Layout(;
        title_text          = "Comparison of Fit Results",
        xaxis_title_text    = "z<sub>Real</sub> / Ω",
        xaxis_dtick         = _dtick,
        yaxis_title_text    = "z<sub>Imag</sub> / Ω",
        yaxis_dtick         = _dtick,
        # --- Fixed Ratio Axes Configuration
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
        # annotations         = _annotations,
    )
    return PlotlyJS.Plot([_trace_data, _trace_Gamry, _trace_Init, _trace_Optim], plt_layout)
end

# --- initial parameter values via "EquivalentCircuits": -------------------------------------------------------------------
q_init, eq_params_init, Z_init =  MyLibOptimizeEquivalentCircuit(circ_strg_ref, frequ_data, Z_data; _max_iter = 10)
# --- reference fit / GamryFit: --------------------------------------------------------------------------------------------
Z_Gamry,   q_Gamry   = simulate_impedance(circ_strg_ref, params_Gamry,   frequ_data, Z_data)
# --- call optimization function: ------------------------------------------------------------------------------------------
if false
    q_opt, eq_params_opt, Z_opt = fit_three_RP_EC(circ_strg_ref, eq_params_init, frequ_data, Z_data, true);
else
    q_opt, eq_params_opt, Z_opt = fit_three_RP_EC(circ_strg_ref, eq_params_init, frequ_data, Z_Gamry, true);
end

# --- plot comparison: -----------------------------------------------------------------------------------------------------
# --- comparizon of qualities: ---------------------------------------------------------------------------------------------
println(@sprintf("q_Gamry: %10.3g, \tq_init: %10.3g, \tq_opt: %10.3g", q_Gamry, q_init, q_opt))

plot_Nyquist(Z_data, Z_Gamry, Z_init, Z_opt, frequ_data)