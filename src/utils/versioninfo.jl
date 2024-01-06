# This file is a part of Groebner.jl. License is GNU GPL v2.

# Adapted from Nemo.jl
# The license is Simplified "2-clause" BSD License
#   https://github.com/Nemocas/Nemo.jl/blob/master/LICENSE.md

import Pkg
import LibGit2
import Dates

const Groebner_uuid = Base.UUID("0b43b601-686d-58a3-8a1c-6623616c7cd4")

deps = Pkg.dependencies()
if !haskey(deps, Groebner_uuid)
    version() = "build"
else
    ver = deps[Groebner_uuid]
    if occursin("/dev/", ver.source)
        version() = VersionNumber("$(ver.version)-dev")
    else
        version() = VersionNumber("$(ver.version)")
    end
end

function versioninfo()
    print("Groebner.jl version $(version())\n")
    groebnerpath = dirname(dirname((@__DIR__)))

    repo = LibGit2.GitRepo(groebnerpath)

    print("Commit ")
    println(string(LibGit2.head_oid(repo))[1:8])
    print("  Author: ")
    commit = LibGit2.GitCommit(repo, LibGit2.head_oid(repo))
    println("$(LibGit2.author(commit).name) <$(LibGit2.author(commit).email)>")
    print("  Date:   ")
    print(Dates.unix2datetime(LibGit2.author(commit).time))

    finalize(repo)

    println()
    println("Switches:")
    println("  invariants_enabled           = $(invariants_enabled())")
    println("  logging_enabled              = $(logging_enabled())")
    println("  performance_counters_enabled = $(performance_counters_enabled())")

    println("Environment:")
    println("  GROEBNER_NO_THREADED = $(get(ENV, "GROEBNER_NO_THREADED", "0"))")

    nothing
end
