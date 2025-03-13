su_include_dir = include_directories('include')

subdir('lib')
subdir('main')

su_lib = library(
    'su_lib',
    su_lib_srcs + su_main_srcs,
    include_directories: su_include_dir,
    dependencies: [cwp_dep, par_dep],
)

su_dep = declare_dependency(
    include_directories: su_include_dir,
    link_with: su_lib,
    dependencies: [cwp_dep, par_dep],
)