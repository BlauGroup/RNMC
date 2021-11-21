{
  description = "High performance Monte Carlo simulator";

  inputs.nixpkgs.url = github:NixOS/nixpkgs/nixos-21.05;

  outputs = { self, nixpkgs }: {

    defaultPackage.x86_64-linux =
      # Notice the reference to nixpkgs here.
      with import nixpkgs { system = "x86_64-linux"; };
      stdenv.mkDerivation {
        name = "RNMC";
        src = self;

        buildInputs = [
          gcc
          clang
          gsl
          sqlite
          sqlitebrowser
          gdb
          valgrind
        ];


        buildPhase = "CC=clang++ ./build.sh";
        installPhase = "mkdir -p $out/bin; mv ./build/* $out/bin";
      };

  };
}
