using BinDeps

@BinDeps.setup

if (@osx? (Base.blas_vendor() == :openblas64) : false)
    aliases = ["libscsdir64"]
else
    aliases = ["libscsdir"]
end

scs = library_dependency("scs", aliases=aliases)

@osx_only begin
    using Homebrew
    provides(Homebrew.HB, "scs", scs, os = :Darwin)
end

version = "1.0.7"

provides(Sources, URI("https://github.com/cvxgrp/scs/archive/v$version.tar.gz"),
    [scs], os=:Unix, unpacked_dir="scs-$version")

# Windows binaries built in Cygwin as follows:
# CFLAGS="-DDLONG -DLAPACK_LIB_FOUND -DBLAS64 -DBLASSUFFIX=_64_" LDFLAGS="-L$HOME/julia/usr/bin -lopenblas" make CC=x86_64-w64-mingw32-gcc out/libscsdir.dll
# mv out bin64
# make clean
# CFLAGS="-DDLONG -DLAPACK_LIB_FOUND" LDFLAGS="-L$HOME/julia32/usr/bin -lopenblas" make CC=i686-w64-mingw32-gcc out/libscsdir.dll
# mv out bin32
provides(Binaries, URI("http://sourceforge.net/projects/juliadeps-win/files/scs-$version.7z"),
    [scs], unpacked_dir="bin$WORD_SIZE", os = :Windows,
    SHA="3bd51b934c1b7bdaa10cf2ceebf2b09fee25b58866072299892d6586f8ced292")

prefix = joinpath(BinDeps.depsdir(scs), "usr")
srcdir = joinpath(BinDeps.depsdir(scs), "src", "scs-$version/")

if VERSION < v"0.4.0-dev+3949"
    libname = "libscsdir.$(Sys.dlext)"
else
    libname = "libscsdir.$(Libdl.dlext)"
end

ldflags = ""
@osx_only begin
    ldflags = "$ldflags -undefined suppress -flat_namespace"
end
cflags = "-DDLONG -DLAPACK_LIB_FOUND"
if Base.blas_vendor() == :openblas64
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

@BinDeps.install [:scs => :scs]
