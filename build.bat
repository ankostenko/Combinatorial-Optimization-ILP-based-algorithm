@echo off
@REM -fsanitize=address -g3 -fno-omit-frame-pointer -fno-optimize-sibling-calls

set src=./src/source.cpp ./src/utils.cpp ./src/config_param_parsing.cpp ./src/win_console_colors.cpp ./src/cycle_extraction.cpp ./src/ilp.cpp ./src/print/printing.cpp ./src/graph.cpp
set output=./bin/ham.exe
set arguments=-Wall -pedantic -Wextra -O0 -g0
set libraries=-lkernel32 -L./external/ -llibscip -lsoplex

set sanitize=%1
if %sanitize%==true set add_args=-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls

clang %src% -o %output% %arguments% --include-directory src/include --include-directory src/include/scip %add_args% %libraries% -std=c++20