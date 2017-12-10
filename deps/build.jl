using BinDeps

@BinDeps.setup

blasvendor = Base.BLAS.vendor()

direct = library_dependency("libscsdir", aliases=["libscsdir64"])
indirect = library_dependency("libscsindir", aliases=["libscsindir64"])

# TODO: Provide both libs in the "scs" Homebrew package.
# if is_apple()
#     using Homebrew
#     provides(Homebrew.HB, "scs", [direct, indirect], os = :Darwin)
# end

version = "2.0.2"
win_version = "2.0.2" # The windows binaries are not consistent with this version yet.

provides(Sources, URI("https://github.com/cvxgrp/scs/archive/v$version.tar.gz"),
    [direct, indirect], os=:Unix, unpacked_dir="scs-$version")

# Windows binaries built in Cygwin as follows:
# CFLAGS="-DDLONG -DCOPYAMATRIX -DUSE_LAPACK -DCTRLC=1 -DBLAS64 -DBLASSUFFIX=_64_" LDFLAGS="-L$HOME/julia/usr/bin -lopenblas64_" make CC=x86_64-w64-mingw32-gcc out/libscsdir.dll
# mv out bin64
# make clean
# CFLAGS="-DDLONG -DCOPYAMATRIX -DUSE_LAPACK -DCTRLC=1" LDFLAGS="-L$HOME/julia32/usr/bin -lopenblas" make CC=i686-w64-mingw32-gcc out/libscsdir.dll
# mv out bin32

# TODO: Provide binaries for windows
# provides(Binaries, URI("https://cache.julialang.org/https://bintray.com/artifact/download/tkelman/generic/scs-$win_version-r2.7z"),
#     [direct, indirect], unpacked_dir="bin$(Sys.WORD_SIZE)", os = :Windows,
#     SHA="62bb4feeb7d2cd3db595f05b86a20fc93cfdef23311e2e898e18168189072d02")

prefix = joinpath(BinDeps.depsdir(direct), "usr")
srcdir = joinpath(BinDeps.depsdir(direct), "src", "scs-$version/")

ldflags = ""
if is_apple()
    ldflags = "$ldflags -undefined suppress -flat_namespace"
end
cflags = "-DCOPYAMATRIX -DDLONG -DUSE_LAPACK -DCTRLC=1"
if blasvendor == :openblas64
    cflags = "$cflags -DBLAS64 -DBLASSUFFIX=_64_"
end
if blasvendor == :mkl
    if Base.USE_BLAS64
        cflags = "$cflags -DMKL_ILP64 -DBLAS64"
        ldflags = "$ldflags -lmkl_intel_ilp64"
    else
        ldflags = "$ldflags -lmkl_intel"
    end
    cflags = "$cflags -fopenmp"
    ldflags = "$ldflags -lmkl_gnu_thread -lmkl_rt -lmkl_core -lgomp"
end

ENV2 = copy(ENV)
ENV2["LDFLAGS"] = ldflags
ENV2["CFLAGS"] = cflags

provides(SimpleBuild,
    (@build_steps begin
        GetSources(direct)
        CreateDirectory(joinpath(prefix, "lib"))

        FileRule(joinpath(prefix, "lib", "libscsdir.$(Libdl.dlext)"), @build_steps begin
            ChangeDirectory(srcdir)
            setenv(`make BLASLDFLAGS= out/libscsdir.$(Libdl.dlext)`, ENV2)
            `mv out/libscsdir.$(Libdl.dlext) $prefix/lib`
        end)

        FileRule(joinpath(prefix, "lib", "libscsindir.$(Libdl.dlext)"), @build_steps begin
            ChangeDirectory(srcdir)
            setenv(`make BLASLDFLAGS= out/libscsindir.$(Libdl.dlext)`, ENV2)
            `mv out/libscsindir.$(Libdl.dlext) $prefix/lib`
        end)

    end), [direct, indirect], os=:Unix)

@BinDeps.install Dict(:libscsdir => :direct, :libscsindir => :indirect)
