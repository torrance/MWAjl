using Revise
using Pkg
using Pkg.Artifacts
using Pkg.BinaryPlatforms
using SHA

Revise.track(Pkg)
Revise.track(Pkg.Artifacts)

artifact_toml = joinpath(@__DIR__, "../Artifacts.toml")

artifacts = [
    ("libbeam.so", "https://github.com/torrance/MWAjl/raw/artifacts/x86_64/libbeam.so.tar.gz"),
    ("libcasacorejl.so", "https://github.com/torrance/MWAjl/raw/artifacts/x86_64/libcasacorejl.so.tar.gz")
]


for (name, url) in artifacts
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
        platform=Linux(:x86_64)
    )
end
