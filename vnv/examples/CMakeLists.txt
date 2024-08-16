# Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
# HYPRE Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

project(HypreExamples)
cmake_minimum_required(VERSION 3.21)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

find_package(HYPRE REQUIRED)
find_package(Injection REQUIRED)
add_executable(ex1 ex1.c)
target_link_libraries(ex1 PUBLIC HYPRE::HYPRE)
link_vnv_executable(ex1 HYPRE_EX1 c)

