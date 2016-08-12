using BinDeps
using Compat

@BinDeps.setup

if VERSION >= v"0.5.0-dev+4679"
    blasvendor = Base.BLAS.vendor()
else
    blasvendor = Base.blas_vendor()
end


if (is_apple() ? (blasvendor == :openblas64) : false)
    aliases = ["libscsdir64"]
else
    aliases = ["libscsdir"]
end

scs = library_dependency("scs", aliases=aliases)

if is_apple()
    using Homebrew
    provides(Homebrew.HB, "scs", scs, os = :Darwin)
end

version = "1.1.8"
win_version = "1.1.5"

provides(Sources, URI("https://github.com/cvxgrp/scs/archive/v$version.tar.gz"),
    [scs], os=:Unix, unpacked_dir="scs-$version")

# Windows binaries built in Cygwin as follows:
# CFLAGS="-DDLONG -DCOPYAMATRIX -DLAPACK_LIB_FOUND -DCTRLC=1 -DBLAS64 -DBLASSUFFIX=_64_" LDFLAGS="-L$HOME/julia/usr/bin -lopenblas64_" make CC=x86_64-w64-mingw32-gcc out/libscsdir.dll
# mv out bin64
# make clean
# CFLAGS="-DDLONG -DCOPYAMATRIX -DLAPACK_LIB_FOUND -DCTRLC=1" LDFLAGS="-L$HOME/julia32/usr/bin -lopenblas" make CC=i686-w64-mingw32-gcc out/libscsdir.dll
# mv out bin32
provides(Binaries, URI("https://cache.julialang.org/https://bintray.com/artifact/download/tkelman/generic/scs-$win_version-r2.7z"),
    [scs], unpacked_dir="bin$(Sys.WORD_SIZE)", os = :Windows,
    SHA="62bb4feeb7d2cd3db595f05b86a20fc93cfdef23311e2e898e18168189072d02")

prefix = joinpath(BinDeps.depsdir(scs), "usr")
srcdir = joinpath(BinDeps.depsdir(scs), "src", "scs-$version/")

libname = "libscsdir.$(Libdl.dlext)"

ldflags = ""
if is_apple()
    ldflags = "$ldflags -undefined suppress -flat_namespace"
end
cflags = "-DCOPYAMATRIX -DDLONG -DLAPACK_LIB_FOUND -DCTRLC=1"
if blasvendor == :openblas64
    cflags = "$cflags -DBLAS64 -DBLASSUFFIX=_64_"
end

ENV2 = copy(ENV)
ENV2["LDFLAGS"] = ldflags
ENV2["CFLAGS"] = cflags

provides(SimpleBuild,
    (@build_steps begin
        GetSources(scs)
        CreateDirectory(joinpath(prefix, "lib"))
        FileRule(joinpath(prefix, "lib", libname), @build_steps begin
            ChangeDirectory(srcdir)
            setenv(`make out/$libname`, ENV2)
            `mv out/$libname $prefix/lib`
        end)
    end), [scs], os=:Unix)

@BinDeps.install Dict([(:scs, :scs)])
