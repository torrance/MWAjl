function parse_model(f::IOStream)::Array{Source}
    parse_version(f)

    sources = Source[]
    while true
        token = gettoken(f)

        if token == "source"
            push!(sources, parse_source(f))
        elseif token == ""
            break
        else
            error("Expected source definition, got: $(token)")
        end
    end
   
    sources
end


function parse_version(f::IOStream)
    # Read version header
    for expected in ["skymodel", "fileformat", "1.1"]
        token = gettoken(f)
        if token != expected
            error("Expected skymodel version 1.1")
        end
    end
end


function parse_source(f::IOStream)::Source
    token = gettoken(f)
    if token != "{"
        error("Expected {, got: $(token)")
    end

    local name::String
    comps = Component[]
    while true
        token = gettoken(f)
        if token == "name"
            name = gettoken(f)
        elseif token == "component"
            push!(comps, parse_component(f))
        elseif token == "}"
            break
        elseif token == ""
            error("Unexpected EOF")
        else
            error("Expected source definition, got $(token)")
        end
    end

    try
        return Source(name, comps)
    catch e
        if isa(e, UndefVarError)
            error("Incomplete source definition")
        else
            rethrow()
        end
    end
end


function parse_component(f::IOStream)::Component
    token = gettoken(f)
    if token != "{"
        error("Expected {, got $(token)")
    end
    
    local position::Position
    local type::String
    local sed::SED
    measurements = Measurement[]
    while true
        token = gettoken(f)
        if token == "position"
            ra = hms2rad(gettoken(f))
            dec = dms2rad(gettoken(f))
            position = Position(ra, dec)
        # TODO: Implement Gaussian type
        elseif token == "type"
            type = gettoken(f)
            if type ∉ ["point"]
                error("Expected component type 'point', got: $(type)")
            end
        elseif token == "sed"
            sed = parse_sed(f)
        elseif token == "measurement"
            push!(measurements, parse_measurement(f))
        elseif token == "}"
            break
        elseif token == ""
            error("Unexpected EOF")
        else
            error("Expected component definition, got: $(token)")
        end
    end

    try
        if length(measurements) > 0
            sort!(measurements, by=x -> x.ν)
            return Component(position, type, Measurements(measurements))
        else
            return Component(position, type, sed)
        end
    catch e
        if isa(e, UndefVarError)
            error("Incomplete component definition")
        else
            rethrow()
        end
    end
end


function parse_measurement(f::IOStream)::Measurement
    token = gettoken(f)
    if token != "{"
        error("Expected {, got $(token)")
    end

    local frequency::Float64
    local I::Float64, Q::Float64, U::Float64, V::Float64
    while true
        token = gettoken(f)
        if token == "frequency"
            try
                frequency = parse(Float64, gettoken(f)) * 1E6
            catch e
                error("Unable to parse frequency value")
            end

            suffix = gettoken(f)
            if suffix != "mhz"
                error("Expected unit MHz, got: $(suffix)")
            end
        elseif token == "fluxdensity"
            suffix = gettoken(f)
            if suffix != "jy"
                error("Expected unit Jy, got: $(suffix)")
            end

            try
                I = parse(Float64, gettoken(f))
                Q = parse(Float64, gettoken(f))
                U = parse(Float64, gettoken(f))
                V = parse(Float64, gettoken(f))
            catch e
                error("Unable to parse flux value")
            end
        elseif token == "}"
            break
        elseif token == ""
            error("Unexpected EOF")
        else
            error("Expected measurement definition, got $(token)")
        end
    end
    
    try
        return Measurement(frequency, [I, Q, U, V])
    catch e
        if isa(e, UndefVarError)
            error("Incomplete measurement definition")
        end
    end
end


function parse_sed(f::IOStream)::SED
    token = gettoken(f)
    if token != "{"
        error("Expected {, got $(token)")
    end

    local frequency::Float64
    local I::Float64, Q::Float64, U::Float64, V::Float64
    local coeffs::Array{Float64}
    while true
        token = gettoken(f)
        if token == "frequency"
            try
                frequency = parse(Float64, gettoken(f)) * 1E6
            catch e
                error("Unable to parse frequency value")
            end

            suffix = gettoken(f)
            if suffix != "mhz"
                error("Expected unit MHz, got: $(suffix)")
            end
        elseif token == "fluxdensity"
            suffix = gettoken(f)
            if suffix != "jy"
                error("Expected unit Jy, got: $(suffix)")
            end

            try
                I = parse(Float64, gettoken(f))
                Q = parse(Float64, gettoken(f))
                U = parse(Float64, gettoken(f))
                V = parse(Float64, gettoken(f))
            catch e
                error("Unable to parse flux value")
            end
        elseif token == "spectral-index"
            coeffs = parse_spectral_index(f)
        elseif token == "}"
            break
        elseif token == ""
            error("Unexpected EOF")
        else
            error("Expected measurement definition, got $(token)")
        end
    end

    try
        return SED(frequency, [I, Q, U, V], coeffs)
    catch e
        if isa(e, UndefVarError)
            error("Incomplete SED definition")
        else
            rethrow
        end
    end
end
        

function parse_spectral_index(f::IOStream)::Array{Float64}
    token = gettoken(f)
    if token != "{"
        error("Expected {, got $(token)")
    end

    coeffs = Float64[]
    while true
        token = gettoken(f)
        if token == "}"
            break
        else
            try
                push!(coeffs, parse(Float64, token))
            catch e
                if isa(e, ArgumentError)
                    error("Unable to parse spectral index value")
                else
                    rethrow()
                end
            end
        end
    end

    if length(coeffs) == 0
        error("At least one spectral index value must be supplied")
    end
    coeffs
end


function hms2rad(hms::String)::Number
    try
        hours, ms = split(hms, 'h', limit=2)
        minutes, s = split(ms, 'm', limit=2)
        seconds, rest = split(s, 's', limit=2)
        
        if rest != ""
            error("rest is not empty")
        end

        hours = Base.parse(Float64, hours)
        minutes = Base.parse(Float64, minutes)
        seconds = Base.parse(Float64, seconds)

        return sign(hours) * deg2rad(abs(hours) * 15 + (minutes / 60) * 15 + (seconds / 3600) * 15)
    catch e
        error("Expected RA of form xxhxxmxxs, got: $(hms)")
    end
end


function dms2rad(dms::String)::Number
    try
        degrees, ms = split(dms, 'd', limit=2)
        minutes, s = split(ms, 'm', limit=2)
        seconds, rest = split(s, 's', limit=2)
        
        if rest != ""
            error("rest is not empty")
        end

        degrees = Base.parse(Float64, degrees)
        minutes = Base.parse(Float64, minutes)
        seconds = Base.parse(Float64, seconds)

        return sign(degrees) * deg2rad(abs(degrees) + (minutes / 60) + (seconds / 3600))
    catch e
        error("Expected Dec of form xxdxxmxxs, got: $(dms)")
    end
end


function gettoken(f::IOStream)::String
    str = ""
    skipchars(isspace, f, linecomment='#')
    while !eof(f)
        char = read(f, Char)
        # First two cases handle quoted strings using " or '.
        # We read the input until the matching quote appears.
        # This does NOT handle any kind of quote escaping.
        if length(str) == 0 && char == '"'
            str = readuntil(f, '"')
            break
        elseif length(str) == 0 && char == '\''
            str = readuntil(f, '\'')
            break
        # Both { and } are full and complete tokens. If we already
        # have a partial token, these indicate the end of this token. Otherwise,
        # just return the brace.
        elseif length(str) > 0 && (char == '{' || char == '}')
            skip(str, -1)  # Rewind the buffer
            break
        elseif char == '{' || char == '}'
            str = string(char)
            break
        # Otherwise, a space indicates the end of a token
        elseif isspace(char)
            break
        else
            str *= char
        end
    end

    lowercase(str)
end