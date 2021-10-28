with (import <nixpkgs> {});

mkShell rec {
  buildInputs = [ gcc
                  clang
                  gsl
                  sqlite
                  sqlitebrowser
                  gdb
                  valgrind
                ];


  # set CPATH for running emacs
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

}
