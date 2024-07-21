{
  description = "High performance Monte Carlo simulator";

  inputs.nixpkgs.url = github:NixOS/nixpkgs/nixos-21.11;

  outputs = { self, nixpkgs }:

    # RNMC is so simple that the build derivation looks exactly the same on
    # all platforms.
    let genericDefaultPackage = systemString:
          with import nixpkgs { system = systemString; };
          stdenv.mkDerivation {
            name = "RNMC";
            src = self;

            buildInputs = [
              clang
              gsl
              sqlite
            ];


            buildPhase = "CC=clang++ ./build.sh";
            installPhase = "mkdir -p $out/bin; mv ./build/* $out/bin";
            doCheck = true;
            checkPhase = "../tests/test.sh";

          };
      in {
        devShell.x86_64-linux =
          with import nixpkgs { system = "x86_64-linux"; };
          mkShell {

            buildInputs = [
              gcc
              clang
              gsl
              (sqlite.override { interactive = true; })
              sqlitebrowser
              gdb
              valgrind
              bintools-unwrapped # gprof 
            ];

            # environment for CLANGD
            # CPATH=$CLANGD_PATH emacs
            CLANGD_PATH = builtins.concatStringsSep ":" [

              # C++ stdlib headers
              "${gcc-unwrapped}/include/c++/10.3.0"
              "${gcc-unwrapped}/include/c++/10.3.0/x86_64-unknown-linux-gnu"

              # libc headers
              "${glibc.dev}/include"

              # compiler specific headers
              "${clang}/resource-root/include"

              "${sqlite.dev}/include"
              "${gsl}/include"
            ];
          };

        defaultPackage = {
          x86_64-linux = genericDefaultPackage "x86_64-linux";
          x86_64-darwin = genericDefaultPackage "x86_64-darwin";
        };

      checks = {
        x86_64-linux.tests = genericDefaultPackage "x86_64-linux";
        x86_64-darwin.tests = genericDefaultPackage "x86_64-darwin";

      };

      };
}
