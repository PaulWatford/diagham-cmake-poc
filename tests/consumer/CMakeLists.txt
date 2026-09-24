# Minimal downstream project: proves an installed DiagHam is usable through
# find_package(DiagHam) alone. Configured, built and run by the
# install.find_package test (tests/InstallConsumer.cmake).
cmake_minimum_required(VERSION 3.16)
project(DiagHamConsumer LANGUAGES CXX)

find_package(DiagHam REQUIRED)

add_executable(diagham_consumer main.cc)
target_link_libraries(diagham_consumer PRIVATE DiagHam::DiagHam)
