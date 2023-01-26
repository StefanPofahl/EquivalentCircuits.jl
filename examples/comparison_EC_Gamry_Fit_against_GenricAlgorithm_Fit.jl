# --- comparizon Gamry Impedance Fit Result and Julia Fit-Result: ----------------------------------------------------------
using RobustModels, PyCall, Printf, PlotlyJS
using StefanPHyCentaImpedance # MylibExpECircJLStrToImpPy
using StefanPHyCentaPlotTools # MyLibPlotlyLayoutHalf!

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
# --- Equivalent Circuit Model: "R1-L2-[P3,R4]-[P5,R6]-[P7,R8]"
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


# --- results of GeneticAlgorithmFit (Version Dec.2022 & Jan.2023): --------------------------------------------------------
# --- Version Dec.2022:
params_GenAlg1 = (R1 = 1.8427875572829713e-7, L2 = 3.842125874199544e-8, 
    P3w = 122.36069899203851,       P3n = 0.01643920897842392, 
    R4  = 4.15262208372949e9, 
    P5w = 3.199017229518061,        P5n = 0.6614015842109463, 
    R6  = 0.010663896457683125, 
    P7w = 1.1381807924643682e10,    P7n = 0.000283697633465696, 
    R8  = 2.927717582694503e8)

# --- Version Jan.2023:
params_GenAlg2 = (R1 = 1.9255835476502468e-10, L2 = 3.842187844027932e-8, 
    P3w = 122.36008204850754, P3n = 0.016437462912654768, 
    R4 = 2.6807001647440516318127683e9, 
    P5w = 3.199108007313649, P5n = 0.6613933758340284, 
    R6 = 0.010664062663627448, 
    P7w = 1.0054456513854383e10, P7n = 0.4969485793197776, 
    R8 = 4.055339967015674e8)


# --- functions ------------------------------------------------------------------------------------------------------------
function simulate_impedance(_circ_strg_ref::String, _params_Gamry::NamedTuple,  
    _frequ_data::Vector{Float64}, _Z_data::Vector{ComplexF64}, _frequ_simu::Vector{Float64})
    # ---
    _circuits = PyCall.pyimport("impedance.models.circuits")
    _circ_str_ImpPy = MylibExpECircJLStrToImpPy(_circ_strg_ref)
    _circuit_ImpPy  = _circuits.CustomCircuit(initial_guess= collect(_params_Gamry), circuit= _circ_str_ImpPy)
    # --- simulate impedance according to given parameters (no optimization):
    _Z_data_simu    = _circuit_ImpPy.predict(_frequ_data, use_initial = true)
    _q_ref          = RobustModels.mean(abs.(_Z_data - _Z_data_simu))
    _Z_simu         = _circuit_ImpPy.predict(_frequ_simu, use_initial = true)
    return _Z_simu, _q_ref
end

function plot_Nyquist(_Z_data::Vector{ComplexF64}, _Z_Gamry::Vector{ComplexF64}, _Z_GenAlg1::Vector{ComplexF64}, _Z_GenAlg2::Vector{ComplexF64}, _frequ_data::Vector{Float64})
    _dtick = 2.5e-3
    s_pts_info = []
    for i_ndx in eachindex(_frequ_data)
        push!(s_pts_info, @sprintf("#:%i, f=%.2fHz", i_ndx, _frequ_data[i_ndx]))
    end
    # --- traces:
    _trace_data     = PlotlyJS.scatter(; x= real(_Z_data),    y=  imag(_Z_data),    name = "measured", text = s_pts_info, mode = "markers")
    _trace_Gamry    = PlotlyJS.scatter(; x= real(_Z_Gamry),   y=  imag(_Z_Gamry),   name = "Gamry")
    _trace_GenAlg1  = PlotlyJS.scatter(; x= real(_Z_GenAlg1), y=  imag(_Z_GenAlg1), name = "GenAlg1")
    _trace_GenAlg2  = PlotlyJS.scatter(; x= real(_Z_GenAlg2), y=  imag(_Z_GenAlg2), name = "GenAlg2")
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
    return PlotlyJS.Plot([_trace_data, _trace_Gamry, _trace_GenAlg1, _trace_GenAlg2], plt_layout)
end

# --- main -----------------------------------------------------------------------------------------------------------------
Z_Gamry,   q_Gamry   = simulate_impedance(circ_strg_ref, params_Gamry,   frequ_data, Z_data, frequ_simu)
Z_GenAlgA, q_GenAlgA = simulate_impedance(circ_strg_ref, params_GenAlg1, frequ_data, Z_data, frequ_simu)
Z_GenAlgB, q_GenAlgB = simulate_impedance(circ_strg_ref, params_GenAlg2, frequ_data, Z_data, frequ_simu)

# --- comparizon of qualities: ---------------------------------------------------------------------------------------------
println(@sprintf("q_Gamry: %10.3g, \tq_GenAlgA: %10.3g, \tq_GenAlgB: %10.3g", q_Gamry, q_GenAlgA, q_GenAlgB))

# --- call plot functions: -------------------------------------------------------------------------------------------------
begin
    fn_plt = raw"c:\tmp\plt\Nyquist_of_fitted_Z.svg"
    hdl_plt = plot_Nyquist(Z_data, Z_Gamry, Z_GenAlgA, Z_GenAlgB, frequ_data); hdl_svg = hdl_plt
    MyLibPlotlyLayoutHalf!(hdl_svg, fn_plt)
    PlotlyJS.savefig(hdl_svg, fn_plt)
end



    