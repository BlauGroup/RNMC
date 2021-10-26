with (import <nixpkgs> {});

mkShell rec {
  buildInputs = [ gcc
                  gsl
                  sqlite
                  gdb
                  valgrind
                  clang
                ];

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

}
