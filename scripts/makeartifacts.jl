using Pkg.Artifacts
using Pkg.BinaryPlatforms
using SHA

artifact_toml = joinpath(@__DIR__, "../Artifacts.toml")

# [ name, platform, lazy, url]
artifacts = [
    ("libbeam.so", Linux(:x86_64), false, "https://github.com/torrance/MWAjl/raw/artifacts/x86_64/libbeam.so.tar.gz"),
    ("libcasacorejl.so", Linux(:x86_64), false, "https://github.com/torrance/MWAjl/raw/artifacts/x86_64/libcasacorejl.so.tar.gz"),
    ("mwabeam", nothing, false, "https://github.com/torrance/MWAjl/releases/download/assets-v1/mwa_full_embedded_element_pattern.h5.tar.gz"),
    ("testdata", nothing, true, "https://github.com/torrance/MWAjl/releases/download/assets-v1/data.tar.gz"),

]

for (name, platform, lazy, url) in artifacts
    download_info = Tuple{String, String}[]

    hash = create_artifact() do artifact_dir
        path = download(url, joinpath(artifact_dir, name * ".tar.gz"))
        open(path) do io
            push!(download_info, (url, bytes2hex(sha256(io))))
        end
        run(`tar -xzf $path -C $artifact_dir`)
        rm(path)
    end

    bind_artifact!(
        artifact_toml,
        name,
        hash,
        force=true,
        download_info=download_info,
        platform=platform,
        lazy=lazy,
    )
end
