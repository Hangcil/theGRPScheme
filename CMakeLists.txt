cmake_minimum_required(VERSION 3.0.0)
project(GRPBookExamples VERSION 0.1.0 LANGUAGES C CXX)

add_executable(GRPBookExample3_11 example3_11.cpp
                               ./GRPForScalarProblems/scalarGRPScheme.cpp)

add_executable(GRPBookExample3_12 example3_12.cpp
                               ./GRPForScalarProblems/scalarGRPScheme.cpp)



