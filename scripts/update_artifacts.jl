#!/usr/bin/env julia

using TOML

function usage()
    println("Usage:")
    println("  julia scripts/update_artifacts.jl <Artifacts.toml> <r8_url> <r8_sha256> <r8_git_tree_sha1> <r16_url> <r16_sha256> <r16_git_tree_sha1>")
    println("")
    println("Example:")
    println("  julia scripts/update_artifacts.jl Artifacts.toml https://.../r8.tar.gz <sha256> <tree> https://.../r16.tar.gz <sha256> <tree>")
end

function set_artifact!(data::Dict{String,Any}, name::String, url::String, sha256::String, tree::String)
    data[name] = Dict(
        "git-tree-sha1" => tree,
        "download" => [Dict("url" => url, "sha256" => sha256)],
    )
end

function main(args)
    if length(args) != 7
        usage()
        error("Expected 7 arguments, got $(length(args)).")
    end

    artifacts_path, r8_url, r8_sha, r8_tree, r16_url, r16_sha, r16_tree = args

    data = if isfile(artifacts_path) && !isempty(strip(read(artifacts_path, String)))
        TOML.parsefile(artifacts_path)
    else
        Dict{String,Any}()
    end

    set_artifact!(data, "spheroidal_backend_r8", r8_url, r8_sha, r8_tree)
    set_artifact!(data, "spheroidal_backend_r16", r16_url, r16_sha, r16_tree)

    open(artifacts_path, "w") do io
        TOML.print(io, data)
    end

    println("Updated $(artifacts_path) with artifact bindings for r8 and r16.")
end

main(ARGS)
