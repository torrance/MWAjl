"""
    Matrix2x2toArray(arr)

Convert Julia's 2x2 matrix to C-ordered 4-array, paying special attention to column major
ordering of Julia versus row major ordered of C. Only used in tests.
"""
function Matrix2x2toArray4(arr::Array{T, N}) where {T, N}
    s = size(arr)
    if N < 2 || s[1] != 2 || s[2] != 2
        throw(ArgumentError("size(arr) must equal (2, 2, ...)"))
    end

    if N == 2
        return reshape(
            permutedims(arr, (2, 1)), 4
        )
    else
        rest = s[3:end]
        return reshape(
            permutedims(arr, (2, 1, 3:N...)), 4, rest...
        )
    end
end

"""
    Array4toMatrix2x2(arr)

Convert C-ordered 4-array to Julia's 2x2 matrix, paying special attention to column major
ordering of Julia versus row major ordered of C. Only used in tests.
"""
function Array4toMatrix2x2(arr::Array{T, N}) where {T, N}
    s = size(arr)
    if N < 1 || s[1] != 4
        throw(ArgumentError("size(arr) must equal (4, ...)"))
    end

    if N == 1
        return permutedims(
            reshape(arr, 2, 2), (2, 1)
        )
    else
        rest = s[2:end]
        return permutedims(
            reshape(arr, 2, 2, rest...), (2, 1, 3:(N+1)...)
        )
    end
end

# We use our own, hardcoded in-place matrix multiplications below, as these
# are faster since we *know* these are 2x2 matrices, plus we incorporate
# the adjoint conjugate into the equation which avoids allocations.
@inline @inbounds function AxB!(C, A, B)
    C[1] = A[1] * B[1] + A[2] * B[3]
    C[2] = A[1] * B[2] + A[2] * B[4]
    C[3] = A[3] * B[1] + A[4] * B[3]
    C[4] = A[3] * B[2] + A[4] * B[4]
end

@inline @inbounds function AxBH!(C, A, B)
    C[1] = A[1] * conj(B[1]) + A[2] * conj(B[2])
    C[2] = A[1] * conj(B[3]) + A[2] * conj(B[4])
    C[3] = A[3] * conj(B[1]) + A[4] * conj(B[2])
    C[4] = A[3] * conj(B[3]) + A[4] * conj(B[4])
end

@inline @inbounds function AHxB!(C, A, B)
    C[1] = conj(A[1]) * B[1] + conj(A[3]) * B[3]
    C[2] = conj(A[1]) * B[2] + conj(A[3]) * B[4]
    C[3] = conj(A[2]) * B[1] + conj(A[4]) * B[3]
    C[4] = conj(A[2]) * B[2] + conj(A[4]) * B[4]
end

@inline @inbounds function plusAxB!(C, A, B)
    C[1] += A[1] * B[1] + A[2] * B[3]
    C[2] += A[1] * B[2] + A[2] * B[4]
    C[3] += A[3] * B[1] + A[4] * B[3]
    C[4] += A[3] * B[2] + A[4] * B[4]
end

@inline @inbounds function plusAxBH!(C, A, B)
    C[1] += A[1] * conj(B[1]) + A[2] * conj(B[2])
    C[2] += A[1] * conj(B[3]) + A[2] * conj(B[4])
    C[3] += A[3] * conj(B[1]) + A[4] * conj(B[2])
    C[4] += A[3] * conj(B[3]) + A[4] * conj(B[4])
end

@inline @inbounds function plusAHxB!(C, A, B)
    C[1] += conj(A[1]) * B[1] + conj(A[3]) * B[3]
    C[2] += conj(A[1]) * B[2] + conj(A[3]) * B[4]
    C[3] += conj(A[2]) * B[1] + conj(A[4]) * B[3]
    C[4] += conj(A[2]) * B[2] + conj(A[4]) * B[4]
end

@inline function AdivB!(C, A, B)
    f =  B[1] * B[4] - B[2] * B[3]
    if f == 0
        throw(SingularException(0))
    end

    C[1] = (A[1] * B[4] - A[2] * B[3]) / f
    C[2] = (A[2] * B[1] - A[1] * B[2]) / f
    C[3] = (A[3] * B[4] - A[4] * B[3]) / f
    C[4] = (A[4] * B[1] - A[3] * B[2]) / f
end

@inline function invA!(A)
    f =  A[1] * A[4] - A[2] * A[3]
    if f == 0
        throw(SingularException(0))
    end

    A[1], A[4] = A[4] / f, A[1] / f
    A[2] = -A[2] / f
    A[3] = -A[3] / f
end
