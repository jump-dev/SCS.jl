using BinDeps

@BinDeps.setup

@unix_only begin
  scs = library_dependency("scs", aliases=["libscsdir"])
end

provides(Sources, URI("https://github.com/cvxgrp/scs/archive/master.zip"),
  [scs], os=:Unix, unpacked_dir="scs-master")

prefix = joinpath(BinDeps.depsdir(scs), "usr")
srcdir = joinpath(BinDeps.depsdir(scs), "src", "scs-master/")


provides(SimpleBuild,
  (@build_steps begin
    GetSources(scs)
    CreateDirectory(joinpath(prefix, "lib"))
    FileRule(joinpath(prefix, "lib", "libscsdir.dylib"), @build_steps begin
      ChangeDirectory(srcdir)
      `cat ${BinDeps.depsdir(scs)}/scs.patch` |> `patch scs.mk`
      `make`
      `mv out/libscsdir.dylib $prefix/lib`
    end)
  end), [scs], os=:Darwin)


provides(SimpleBuild,
  (@build_steps begin
    GetSources(scs)
    CreateDirectory(joinpath(prefix, "lib"))
    FileRule(joinpath(prefix, "lib", "libscsdir.so"), @build_steps begin
      ChangeDirectory(srcdir)
      `cat ${BinDeps.depsdir(scs)}/scs.patch` |> `patch scs.mk`
      `make`
      `mv out/libscsdir.so $prefix/lib`
    end)
  end), [scs], os=:Unix)

@BinDeps.install [:scs => :scs]
