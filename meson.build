project(
    'su',
    'c', 'cpp',
    version : '4426',
    default_options: ['cpp_std=c++17'],
    license: 'GPLv3',
)

cc = meson.get_compiler('c')

# String literals are `const char[]` in C++, and the SU code gives them to `char *` parameters everywhere (err("..."), ...).
# GCC and clang allow that; MSVC needs to be told to.
cpp_compiler = meson.get_compiler('cpp')
if cpp_compiler.get_id() == 'msvc'
  add_project_arguments('/Zc:strictStrings-', language : 'cpp')
else
  add_project_arguments(cpp_compiler.get_supported_arguments('-Wno-write-strings', '-Wno-writable-strings'), language : 'cpp')
endif

subdir('src')