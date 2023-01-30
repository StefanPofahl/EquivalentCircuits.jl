# --- ----------------------------------------------------------------------------------------------------------------------
# https://www.pcre.org/current/doc/html/pcre2syntax.html
# https://perldoc.perl.org/perlre#Modifiers
# --- ----------------------------------------------------------------------------------------------------------------------
using CSV, Printf
fn_path, _ = splitdir(@__FILE__())
fn_path = joinpath(fn_path, "data")
fn_data = joinpath(fn_path, "electrolyse_stack.DTA")
if ~isfile(fn_data)
    error("No file found: \"", fn_data, "\"!")
end

# --- use regular expression to find the part of the text of interest:
    # --- ------------------------------------------------------------------------------------------------------------------
    # internet sites about regular expressions:
    # https://www.pcre.org/current/doc/html/pcre2syntax.html
    # https://perldoc.perl.org/perlre#Modifiers
    # --- some examples of Regular expressions: ----------------------------------------------------------------------------
    # \d:   decimal digit ||  \s: white space ||  \D: not a decimal digit  || .  any character except newline 
    # \V:   character that is not a vertical white space character (includes as well punctuation chars)
    # [a-zA-Z] all a-Z Characters.
    # [a-zA-ZäöüÄÖÜ] all a-Z Characters and also äöüÄÖÜ
    # +     one or more, greedy   || ^ Anchor ^a matches only if a is the first char 
    # \\:   Escape special letter, e.g. \\[ to find left squared bracket
    # ----------------------------------------------------------------------------------------------------------------------
# --- search begin of blocks and end of blocks
if isfile(fn_data)
    txt_lines = readlines(fn_data, keep=true)
    n_lines = size(txt_lines)[1]
    idx_header_line = nothing
    char_dicimal_delim = '.'
    s_header = ""
    r_search_pattern_first_line = Regex("^ZCURVE")
    r_comma_in_numbers = Regex("\t\\d+,\\d+")
    idx_line = 0
    for i_line in txt_lines
        global idx_line, idx_header_line, char_dicimal_delim, s_header
        idx_line += 1
        if occursin(r_search_pattern_first_line, i_line)
            idx_header_line = idx_line + 3
            @info("first line found :-) #: ", idx_header_line)
        end
        if idx_line == idx_header_line
            println("1st data Line: ", i_line)
        end
        if ~isnothing(idx_header_line) && (idx_line == idx_header_line + 1)
            if occursin(r_comma_in_numbers, i_line)
                char_dicimal_delim = ','
                @info("Comma is decimal delimiter!")
            end
        end
    end
end

if isnothing(idx_header_line)
    error("Begin of data block not found!")
else
    s_header = ["", "Point", "Time", "Frequ", "Z_real", "Z_imag", "Z_sig", "Z_mod", "Z_phz", "I_dc", "V_dc", "IE_range"]
    filecontend = CSV.File(fn_data; 
        header = s_header,
        skipto = idx_header_line, 
        decimal = char_dicimal_delim, 
        delim = "\t");
end

n_colums = size(filecontend.names)[1]
if n_colums != size(s_header)[1]
    error("Column Missmatch")
end

names_ = filecontend.names
println(@sprintf("%s, ", String.(names_)))

frequ_data = filecontend[:Frequ]
Z_data = ComplexF64.(filecontend[:Z_real], filecontend[:Z_imag])

function _MyLibReadGamryDTA(fn_DTA::AbstractString)
    if ~isfile(fn_DTA)
        error("No file found: \"", fn_DTA, "\"!")
    end
    txt_lines = readlines(fn_data, keep=true)
    idx_header_line = nothing
    char_dicimal_delim = '.'
    r_search_pattern_first_line = Regex("^ZCURVE")
    r_comma_in_numbers = Regex("\t\\d+,\\d+")
    idx_line = 0
    for i_line in txt_lines
        idx_line += 1
        if occursin(r_search_pattern_first_line, i_line)
            idx_header_line = idx_line + 3
            @info("first line found :-) #: ", idx_header_line)
        end
        if idx_line == idx_header_line
            println("1st data Line: ", i_line)
        end
        if ~isnothing(idx_header_line) && (idx_line == idx_header_line + 1)
            if occursin(r_comma_in_numbers, i_line)
                char_dicimal_delim = ','
                @info("Comma is decimal delimiter!")
            end
        end
    end
    # --- extract data
    if isnothing(idx_header_line)
        error("Begin of data block not found!")
    else
        # --- original: Pt          Time	Freq	 Zreal	   Zimag	 Zsig	      Zmod	   Zphz	    Idc	    Vdc	    IERange
        s_header = ["", "Point",    "Time", "Frequ", "Z_real", "Z_imag", "Amplitude", "Z_mod", "Z_phz", "I_dc", "V_dc", "IE_range"]
        filecontend = CSV.File(fn_data; 
            header = s_header,
            skipto = idx_header_line, 
            decimal = char_dicimal_delim, 
            delim = "\t");
    end
    # --- check column missmatch:
    n_colums = size(filecontend.names)[1]
    if n_colums != size(s_header)[1]
        error("Column Missmatch")
    end
    _frequ_data = filecontend[:Frequ]
    _Z_data = ComplexF64.(filecontend[:Z_real], filecontend[:Z_imag])
    idx_sorted = sortperm(_frequ_data)
    # ---
    return _frequ_data[idx_sorted], _Z_data[idx_sorted]
end

frequ_data_all, Z_data_all = _MyLibReadGamryDTA(fn_data)


