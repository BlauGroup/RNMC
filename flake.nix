{
  description = "High performance Monte Carlo simulator";

  inputs.nixpkgs.url = github:NixOS/nixpkgs/nixos-21.05;

  outputs = { self, nixpkgs }: {

    devShell.x86_64-linux =
      with import nixpkgs { system = "x86_64-linux"; };
      mkShell {

        buildInputs = [
          gcc
          clang
          gsl
          sqlite
          sqlitebrowser
          gdb
          valgrind
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

    defaultPackage.x86_64-linux =
      with import nixpkgs { system = "x86_64-linux"; };
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

      };

  };
}
