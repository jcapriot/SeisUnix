par_lib_src = [
    'atopkge.cpp',
     #'docpkge.c',
    'ealloc.cpp',
    'errpkge.cpp',
    # 'filestat.c',
    'fractal.cpp',
    'getpars.cpp',
    #'lincoeff.c',
    'minfunc.cpp',
    'modeling.cpp',
    'refaniso.cpp',
    'rke.cpp',
    'smooth.cpp',
    'subcalls.cpp',
    #'syscalls.c',
    'taup.cpp',
    'upweik.cpp',
    'VND.cpp',
    'wtlib.cpp',
]

par_lib = library(
    'par',
    par_lib_src,
    include_directories: par_include_dir,
    dependencies: su_cwp_dep
)

su_par_dep = declare_dependency(
    include_directories: par_include_dir,
    link_with: par_lib,
    dependencies: su_cwp_dep
)