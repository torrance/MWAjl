using Documenter, MWAjl

makedocs(
    sitename="MWAjl Documentation",
    modules=[MWAjl],
    format=Documenter.HTML(
        assets=[
            joinpath("assets", "styles.css"),
        ],
    ),
)