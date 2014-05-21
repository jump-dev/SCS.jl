using BinDeps

@BinDeps.setup

@unix_only begin
  scs = library_dependency("scs", aliases=["libscsdir"])
end

provides(Sources, URI("https://github.com/cvxgrp/scs/archive/1.0.2.zip"),
  [scs], os=:Unix, unpacked_dir="scs-1.0.2")

prefix = joinpath(BinDeps.depsdir(scs), "usr")
srcdir = joinpath(BinDeps.depsdir(scs), "src", "scs-1.0.2/")


provides(SimpleBuild,
  (@build_steps begin
    GetSources(scs)
    CreateDirectory(joinpath(prefix, "lib"))
    FileRule(joinpath(prefix, "lib", "libscsdir.dylib"), @build_steps begin
      ChangeDirectory(srcdir)
      `cat ${BinDeps.depsdir(scs)}/make-dylib.patch` |> `patch Makefile`
      `cat ${BinDeps.depsdir(scs)}/scs-fpic.patch` |> `patch scs.mk`
      `make libscsdir.dylib`
      `mv libscsdir.dylib $prefix/lib`
    end)
  end), [scs], os=:Darwin)


provides(SimpleBuild,
  (@build_steps begin
    GetSources(scs)
    CreateDirectory(joinpath(prefix, "lib"))
    FileRule(joinpath(prefix, "lib", "libscsdir.so"), @build_steps begin
      ChangeDirectory(srcdir)
      `cat ${BinDeps.depsdir(scs)}/make-so.patch` |> `patch Makefile`
      `cat ${BinDeps.depsdir(scs)}/scs-fpic.patch` |> `patch scs.mk`
      `make libscsdir.so`
      `mv libscsdir.so $prefix/lib`
    end)
  end), [scs], os=:Unix)

@BinDeps.install [:scs => :scs]
