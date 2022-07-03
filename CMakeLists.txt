cmake_minimum_required (VERSION 3.8)

add_library (hexa_audio INTERFACE)

target_include_directories (hexa_audio INTERFACE include/)

target_compile_features(hexa_audio INTERFACE cxx_std_17)
