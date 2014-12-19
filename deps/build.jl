using BinDeps

@BinDeps.setup

@unix_only begin
    scs = library_dependency("scs", aliases=["libscsdir64","libscsdir"])
end

@osx_only begin
    using Homebrew
    provides(Homebrew.HB, "scs", scs, os = :Darwin)
end

version = "1.0.7"

provides(Sources, URI("https://github.com/cvxgrp/scs/archive/v$version.tar.gz"),
    [scs], os=:Unix, unpacked_dir="scs-$version")

prefix = joinpath(BinDeps.depsdir(scs), "usr")
srcdir = joinpath(BinDeps.depsdir(scs), "src", "scs-$version/")

libname = "libscsdir.$(Sys.dlext)"

@osx_only begin
    ldflags = "-undefined suppress -flat_namespace"
    ENV["LDFLAGS"] = ldflags
end
cflags = "-DDLONG -DLAPACK_LIB_FOUND"
if Base.blas_vendor() == :openblas64
    cflags = "$cflags -DBLAS64 -DBLASSUFFIX=_64_"
end
ENV["CFLAGS"] = cflags

provides(SimpleBuild,
    (@build_steps begin
        GetSources(scs)
        CreateDirectory(joinpath(prefix, "lib"))
        FileRule(joinpath(prefix, "lib", libname), @build_steps begin
            ChangeDirectory(srcdir)
            `make out/$libname`
            `mv out/$libname $prefix/lib`
        end)
    end), [scs], os=:Unix)

@BinDeps.install [:scs => :scs]
