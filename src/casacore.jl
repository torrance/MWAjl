using Pkg.Artifacts

mutable struct Table
    ptr::Ptr{Cvoid}
    function Table(ptr::Ptr{Cvoid})
        tbl = new(ptr)
        finalizer(close, tbl)
        return tbl
    end
end

mutable struct Frame
    ptr::Ptr{Cvoid}
    function Frame(ptr::Ptr{Cvoid})
        frame = new(ptr)
        finalizer(del, frame)
        return frame
    end
end

@enum CasaError OK TableNoFile ArraySlicerError TableError

@enum(TypeEnum,
      TpBool, TpChar, TpUChar, TpShort, TpUShort, TpInt, TpUInt,
      TpFloat, TpDouble, TpComplex, TpDComplex, TpString, TpTable,
      TpArrayBool, TpArrayChar, TpArrayUChar, TpArrayShort, TpArrayUShort,
      TpArrayInt, TpArrayUInt, TpArrayFloat, TpArrayDouble, TpArrayComplex,
      TpArrayDComplex, TpArrayString, TpRecord, TpOther, TpQuantity,
      TpArrayQuantity, TpInt64, TpArrayInt64, TpNumberOfTypes)

const type2str = Dict(Bool             => :boolean,  Int32     => :int,
                      Float32          => :float,    Float64   => :double,
                      Complex{Float32} => :complex,  String    => :string)

const enum2type = Dict(TpBool    => Bool,             TpArrayBool    => Array{Bool},
                       TpInt     => Int32,            TpArrayInt     => Array{Int32},
                       TpFloat   => Float32,          TpArrayFloat   => Array{Float32},
                       TpDouble  => Float64,          TpArrayDouble  => Array{Float64},
                       TpComplex => Complex{Float32}, TpArrayComplex => Array{Complex{Float32}},
                       TpString  => String,           TpArrayString  => Array{String})

const libcasacore = joinpath(artifact"libcasacorejl.so", "libcasacorejl.so")

function Table(path::String)::Table
    err = Ref{Cint}(0)
    ptr = ccall((:table_open, libcasacore), Ptr{Cvoid}, (Cstring, Ptr{Cint}), path, err)
    if CasaError(err[]) == TableNoFile
        throw(ArgumentError("Could not open table: $path"))
    end

    return Table(ptr)
end

function close(tbl::Table)
    # This cannot be a ccall, as no task switch is allowed during finalizers.
    # Otherwise, this results in difficult to diagnose concurrency violation errors.
    ccall((:table_close, libcasacore), Cvoid, (Ptr{Cvoid},), tbl.ptr)
end

function taql(query::String, tbl::Table)
    ptr = ccall((:taql, libcasacore), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring), tbl.ptr, query)
    # TODO: Handle errors
    return Table(ptr)
end

function column_info(tbl::Table, colname::String)
    exists = ccall((:column_exists, libcasacore), Bool, (Ptr{Cvoid}, Cstring), tbl.ptr, colname)
    if !exists
        throw(ArgumentError(string("Column ", colname, " does not exist")))
    end

    ndim, typeid = Ref{Cint}(0), Ref{Cint}(0)
    shape_ptr = ccall((:column_info, libcasacore), Ptr{Csize_t}, (Ptr{Cvoid}, Cstring, Ptr{Cint}, Ptr{Cint}), tbl.ptr, colname, typeid, ndim)
    shape = Tuple(convert(Vector{Int}, unsafe_wrap(Vector{Csize_t}, shape_ptr, ndim[], own=true)))
    return enum2type[TypeEnum(typeid[])], shape
end

#=
Example: column(table, "DATA", blc=[1, 24], trc=[4, 24])
will return full 2x2 matrices from channel 24

Note: slices are one indexed.
=#
function column(tbl::Table, colname::String; blc::Vector{Int} = Int[], trc::Vector{Int} = Int[])
    if length(blc) != length(trc)
        throw(ArgumentError("blc and trc must have same length"))
    end
    # blc and trc are one indexed inside Julia, and zero indexed in casacore
    blc .-= 1
    trc .-= 1

    T, _ = column_info(tbl, colname)
    return column(tbl, colname, blc, trc, T)
end

for T in (Bool, Int32, Float32, Float64, Complex{Float32})
    cname = String(Symbol(:get_column_, type2str[T]))
    @eval function column(tbl::Table, colname::String, blc::Vector{Int}, trc::Vector{Int}, ::Type{$T})::Array{$T}
        ndim, shape_ptr, err = Ref{Cint}(0), Ref{Ptr{Csize_t}}(0), Ref{Int}(0)
        data_ptr = ccall(
            ($cname, libcasacore),
            Ptr{$T},
            (Ptr{Cvoid}, Cstring, Ptr{Cint}, Ptr{Ptr{Csize_t}}, Cint, Ptr{Csize_t}, Ptr{Csize_t},Ptr{Int}),
            tbl.ptr, colname, ndim, shape_ptr, length(blc), blc, trc, err
        )
        if CasaError(err[]) == ArraySlicerError
            _, shape = column_info(tbl, colname)
            throw(ArgumentError(string("Slices (", blc, " and  ", trc, ") are invalid for column shape ", shape)))
        elseif CasaError(err[]) == TableError
            error(string("TableError trying to access column ", colname))
        end
        shape = Tuple(unsafe_wrap(Vector{Csize_t}, shape_ptr[], ndim[], own=true))
        return unsafe_wrap(Array{$T}, data_ptr, shape, own=true)
    end
end

function Frame(mjd::Float64, lon::Float64, lat::Float64)
    ptr = ccall((:frame_new, libcasacore), Ptr{Cvoid}, (Cdouble, Cdouble, Cdouble), mjd, lon, lat)
    return Frame(ptr)
end

function del(frame::Frame)
    ccall((:frame_del, libcasacore), Cvoid, (Ptr{Cvoid},), frame.ptr)
end

function radec2altaz(pos::Position, frame::Frame)
    alt = Ref{Cdouble}(0)
    az = Ref{Cdouble}(0)
    ccall((:radec_to_altaz, libcasacore), Cvoid, (Cdouble, Cdouble, Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), pos.ra, pos.dec, frame.ptr, alt, az)
    return alt[], az[]
end