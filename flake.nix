{
  description = "High performance Monte Carlo simulator";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/master;
    flake-compat = {
      url = github:edolstra/flake-compat;
      flake = false;
    };
    mini-compile-commands = {
      url = github:danielbarter/mini_compile_commands;
      flake = false;
    };
  };

  outputs = { self, nixpkgs, flake-compat, mini-compile-commands }:

    # RNMC is so simple that the build derivation looks exactly the same on
    # all platforms.
    let genericDefaultPackage = systemString:
          with import nixpkgs { system = systemString; };
          clang13Stdenv.mkDerivation {
            name = "RNMC";
            src = self;

            buildInputs = [
              gsl
              sqlite
            ];


            buildPhase = "./build.sh";
            installPhase = "mkdir -p $out/bin; mv ./build/* $out/bin";
            doCheck = true;
            checkPhase = "tests/test.sh";
          };

      in {
        devShell.x86_64-linux =
          with import nixpkgs { system = "x86_64-linux"; };
          # (mkShell.override { stdenv = (callPackage mini-compile-commands {}).wrap clangStdenv; }) {
          (mkShell.override { stdenv = clangStdenv; }) {
            buildInputs = [
              gcc
              clang
              gsl
              (sqlite.override { interactive = true; })
              sqlitebrowser
              gdb
              valgrind
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
