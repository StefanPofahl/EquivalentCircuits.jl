
using PyCall
using PlotlyJS
using Printf, RobustModels, EquivalentCircuits
import Conda, Pkg
# --- install master of EquivalentCircuits: --------------------------------------------------------------------------------
# change to Pkg-mode via "]":
# pkg> add https://github.com/MaximeVH/EquivalentCircuits.jl.git
# --------------------------------------------------------------------------------------------------------------------------
# impedance.py:
# impedance.py is a well known python package for equivalent circuit parameter optimization:
# repository: https://github.com/ECSHackWeek/impedance.py
# --------------------------------------------------------------------------------------------------------------------------
# Conda.jl:
# repository: https://github.com/JuliaPy/Conda.jl
# --- install "impedance" in the default main julia-conda environment
# Conda.add("impedance") 
# --- if instalation via "repo.anaconda.com" maintained by "Anaconda, Inc" fails try "Miniforge":
# julia> ENV["CONDA_JL_USE_MINIFORGE"] = "1"
# change to Pkg-mode via "]":
# pkg> build Conda
# --- it might be that it is necessary to specify a "channel":
# add channel: 
# import Conda; Conda.add_channel("pypi")
# import Conda; Conda.update()
# import Conda; Conda.add("impedance"; channel="pypi")
# --------------------------------------------------------------------------------------------------------------------------
# --- install impedance under Linux via CondaPkg (if installation via Conda failed):
# --- remove old installation:
# import Conda; Conda.rm("impedance")
# import Conda; Conda.pip("uninstall", "impedance")
# using CondaPkg; CondaPkg.rm_pip("impedance")
# CondaPkg.rm_channel("pypi")
# --- re-install via CondaPkg:
# using CondaPkg; CondaPkg.add("impedance")
# using CondaPkg; CondaPkg.add_pip("impedance")
# using CondaPkg; CondaPkg.add(impedance; version="1.4.1", channel="pypi")
# --------------------------------------------------------------------------------------------------------------------------
import Pkg; Pkg.status("EquivalentCircuits")
# --- load python-package: "impedance.py" ----------------------------------------------------------------------------------
circuits = PyCall.pyimport("impedance.models.circuits")
# --- alternative import command:
# circuits = PyCall.pyimport_conda("impedance.models.circuits", "pypi")

# --- example data: --------------------------------------------------------------------------------------------------------
frequ_data = [0.0199553, 0.0251206, 0.0316296, 0.0398258, 0.0501337, 0.0631739, 0.0794492, 0.1001603, 0.1260081, 
    0.1588983, 0.2003205, 0.2520161, 0.316723, 0.400641, 0.5040323, 0.6334459, 0.7923428, 0.999041, 1.266892, 
    1.584686, 1.998082, 2.504006, 3.158693, 3.945707, 5.008013, 6.317385, 7.944915, 9.93114, 12.40079, 15.625, 
    19.86229, 24.93351, 31.25, 38.42213, 50.22321, 63.3446, 79.00281, 100.4464, 125.558, 158.3615, 198.6229, 
    252.4038, 315.5048, 397.9953, 505.5147, 627.7902, 796.875, 998.264, 1265.625, 1577.524, 1976.103, 2527.573]
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

# --- local function: -------------------------------------------------------------------------------------------------------
function MylibExpECircJLStrToImpPy(_circuit_str::String)
    left_round_brackets =  replace(_circuit_str, "[" => "(")
    round_brackets      =  replace(left_round_brackets, "]" => ")")
    char_P_replaced     = replace(round_brackets, "P" => "CPE") # P = constant phase element
    return replace(char_P_replaced, "(" => "p(")
end

function _MyLibEqCircNamedTuple(_circuit_str::String, _param_vec::Vector{<:Number})
    _r_search_pat = r"[RLCPW]{1}[0-9]+"
    _vec_param_symbol = Vector{Symbol}(undef, 0)
    idx_start = 1
    for i_ = 1:length(_circuit_str) - 1
        # global idx_start
        match_result    = match(_r_search_pat, _circuit_str, idx_start)
        if isnothing(match_result)
            break
        else
            if String(match_result.match)[1] == 'P'
                push!(_vec_param_symbol, Symbol(string(match_result.match, "w")))
                push!(_vec_param_symbol, Symbol(string(match_result.match, "n")))
            else
                push!(_vec_param_symbol, Symbol(match_result.match))
            end
            _, idx_start = findfirst(match_result.match, _circuit_str)
            idx_start += 1
        end
    end
    # ---
    return NamedTuple.(zip.(Ref(_vec_param_symbol), zip(_param_vec...)))
end

function _MyLibOptimizeEquivalentCircuit(_circ_strg::String, _z_vec::Vector{ComplexF64}, _frequ::Vector{<:Number}; 
    _max_iter::Int=20 )
    _stag_max = max(1, min(_max_iter, trunc(Int, 0.5 * _max_iter)))
    _circfunc = EquivalentCircuits.circuitfunction(_circ_strg)
    _i_Q = +Inf; i_stagnation = 0
    _Q_best = +Inf; _circ_params_best = []; _Z_simulated_best = []
    for i_ = 1:_max_iter
        _circuit_params_EqCirc  = EquivalentCircuits.parameteroptimisation(_circ_strg, _z_vec, _frequ) 
        _Z_simulated_EqCirc     = EquivalentCircuits.simulateimpedance_noiseless(_circfunc, _circuit_params_EqCirc, _frequ)
        _Q_EqCirc               = RobustModels.mean(abs.(_z_vec - _Z_simulated_EqCirc))
        # --- optimize via Impedance.py / "ImpPy" -----------------------------------------------------
        _circ_str_ImpPy  = MylibExpECircJLStrToImpPy(_circ_strg)
        _initial_ImpPy   = collect(_circuit_params_EqCirc)
        _circuit_ImpPy   = circuits.CustomCircuit(initial_guess= _initial_ImpPy, circuit= _circ_str_ImpPy)
        try
            _circuit_ImpPy.fit(_frequ, _z_vec)         
        catch _error_msg
            @warn(string("ImpPy: Circuit Fit Failed"), exception = (_error_msg, catch_backtrace()))
        end
        if isempty(_circuit_ImpPy.parameters_)
            _Q_ImpPy = +Inf
        else
            _circuit_params_ImpPy   = _MyLibEqCircNamedTuple(_circ_strg, _circuit_ImpPy.parameters_)
            _Z_simulated_ImpPy      = EquivalentCircuits.simulateimpedance_noiseless(_circfunc, _circuit_params_ImpPy, _frequ)
            _Q_ImpPy                = RobustModels.mean(abs.(_z_vec - _Z_simulated_ImpPy))
        end
        if min(_Q_ImpPy, _Q_EqCirc) < _Q_best
            if _Q_ImpPy < _Q_EqCirc
                _Q_best = _Q_ImpPy
                _circ_params_best = _circuit_params_ImpPy
                _Z_simulated_best = _Z_simulated_ImpPy
            else
                _Q_best = _Q_EqCirc
                _circ_params_best = _circuit_params_EqCirc
                _Z_simulated_best = _Z_simulated_EqCirc
            end
            i_stagnation = 0
        else
            i_stagnation += 1
        end
        if i_stagnation > _stag_max
            break
        end
        # println("i_: ", i_, "/", _max_iter, ", i_stag: ", i_stagnation, "/", _stag_max, ", \tQ_best: ", _Q_best, 
        # ", \tΔQ=(Q_best-Q): ", _Q_best - min(_Q_ImpPy, _Q_EqCirc), ", ΔQ=(Q_ImpPy-Q_EqCirc): ", _Q_ImpPy - _Q_EqCirc)
        println(@sprintf("i_: %2i/%i, i_stag: %2i/%i, \tQ_best: %9.3g, \tΔQ=(Q_best-Q): %9.3g, ΔQ=(Q_ImpPy-Q_EqCirc): %9.3g", 
        i_, _max_iter, i_stagnation, _stag_max, _Q_best, _Q_best - min(_Q_ImpPy, _Q_EqCirc) , _Q_ImpPy - _Q_EqCirc))
    end
    return _Q_best, _circ_params_best, _Z_simulated_best 
end

# --- plot nyquist: --------------------------------------------------------------------------------------------------------
function plot_nyquist(_title_str::String, _Circuit_str::String, _frequ_data::Vector{<:Number}, _Z_data::Vector{ComplexF64}, 
    _Z_simulated_data::Vector{ComplexF64}, _data_range::Union{StepRange, Nothing}=nothing)
    # ---
    if isnothing(_data_range)
        _data_range = range(1, step = 1, stop = length(_frequ_data))
    end
    s_pts_info = Vector{String}(undef, 0)
    for i_ndx in eachindex(_frequ_data)
        push!(s_pts_info, @sprintf("#:%i, f=%.2fHz", _data_range[i_ndx], _frequ_data[i_ndx]))
    end
    # ---
    trace_measured   = PlotlyJS.scatter(; x= real(_Z_data), y=  imag(_Z_data), name = "measured", text = s_pts_info, mode = "markers")
    trace_simulated  = PlotlyJS.scatter(; x= real(_Z_simulated_data), y= imag(_Z_simulated_data), name = "simulated",
                                    text = s_pts_info, mode = "markers")
    # --- annotations:
    _annotations = [PlotlyJS.attr(;
    xref        ="paper",           yref    = "paper",
    x           = 0.1,              y       = 0.9,
    text        = _Circuit_str,
    xanchor     = "left",
    yanchor     = "top",
    font_family = "PT Sans Narrow, monospace",
    showarrow   = false)]
    # ---
    plt_layout = PlotlyJS.Layout(;
        title_text          = _title_str,
        xaxis_title_text    = "z<sub>Real</sub>",
        xaxis_dtick         = 1000,
        yaxis_title_text    = "z<sub>Imag</sub>",
        yaxis_dtick         = 1000,
        # --- Fixed Ratio Axes Configuration
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
        annotations         = _annotations,
    )
    return PlotlyJS.Plot([trace_measured, trace_simulated], plt_layout)
end

# --- main -----------------------------------------------------------------------------------------------------------------
circ_strg = "R1-L2-[P3,R4]-[P5,R6]-[P7,R8]"
q_, param_EQ_, Z_simulated_ = _MyLibOptimizeEquivalentCircuit(circ_strg, Z_data, frequ_data; _max_iter = 20);
println("optimzed parameters: ")
println(param_EQ_, "\n", '-'^100)
plot_nyquist("Comparison Measured <-> Simulated", string(circ_strg, ", Q: ", q_), frequ_data, Z_data, Z_simulated_)
