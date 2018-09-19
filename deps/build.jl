using BinaryProvider # requires BinaryProvider 0.3.0 or later

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))
products = [
    LibraryProduct(prefix, ["libscsindir"], :indirect),
    LibraryProduct(prefix, ["libscsdir"], :direct),
]

# Download binaries from hosted location
bin_prefix = "https://github.com/JuliaOpt/SCSBuilder/releases/download/v2.0.2"

# Listing of files generated by BinaryBuilder:
download_info = Dict(
    Linux(:aarch64, :glibc) => ("$bin_prefix/SCSBuilder.v2.0.2.aarch64-linux-gnu.tar.gz", "e5ec3846d94f91558d49ab6c29a96da62fd19e2d932be6d84145b6bf8afc4a5a"),
    Linux(:armv7l, :glibc, :eabihf) => ("$bin_prefix/SCSBuilder.v2.0.2.arm-linux-gnueabihf.tar.gz", "c5ee26dd59e4473787323d3839fa4dd298515eb7df93c42490515ceba47bcabf"),
    Linux(:i686, :glibc) => ("$bin_prefix/SCSBuilder.v2.0.2.i686-linux-gnu.tar.gz", "94ea10faa19378de705a2005c9e67d1048b1beeff26b616fd7968119c63a50c1"),
    Windows(:i686) => ("$bin_prefix/SCSBuilder.v2.0.2.i686-w64-mingw32.tar.gz", "dfb70647ef06fca7652d13ef367914531a1d69d269cc7dd73b0c61324ef1f754"),
    Linux(:powerpc64le, :glibc) => ("$bin_prefix/SCSBuilder.v2.0.2.powerpc64le-linux-gnu.tar.gz", "f83b3f1955d4ab0baf19a8ee0d8ec42f3e88e0013eceac8b1560617c20c3aefd"),
    MacOS(:x86_64) => ("$bin_prefix/SCSBuilder.v2.0.2.x86_64-apple-darwin14.tar.gz", "60cb2bd8adf5d5b98122a2b22deb1e3dc33c446018195ee7fecf4b711870d8ea"),
    Linux(:x86_64, :glibc) => ("$bin_prefix/SCSBuilder.v2.0.2.x86_64-linux-gnu.tar.gz", "4ce7e12dde0bd51500b0e648a9267996692bb503f970f608fb48e33974799d68"),
    FreeBSD(:x86_64) => ("$bin_prefix/SCSBuilder.v2.0.2.x86_64-unknown-freebsd11.1.tar.gz", "59a360e7360114310c0640a13cb3a3cdeb25444ae0ee2be8bf15bf92671969fc"),
    Windows(:x86_64) => ("$bin_prefix/SCSBuilder.v2.0.2.x86_64-w64-mingw32.tar.gz", "ab193056a9d326cf518e0f944010020089246db553aef03dfc64be8ba78ecb59"),
)

# Install unsatisfied or updated dependencies:
unsatisfied = any(!satisfied(p; verbose=verbose) for p in products)
if haskey(download_info, platform_key())
    url, tarball_hash = download_info[platform_key()]
    if unsatisfied || !isinstalled(url, tarball_hash; prefix=prefix)
        # Download and install binaries
        install(url, tarball_hash; prefix=prefix, force=true, verbose=verbose)
    end
elseif unsatisfied
    # If we don't have a BinaryProvider-compatible .tar.gz to download, complain.
    # Alternatively, you could attempt to install from a separate provider,
    # build from source or something even more ambitious here.
    error("Your platform $(triplet(platform_key())) is not supported by this package!")
end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)
