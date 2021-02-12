@echo off

set src=./src/test.cpp ./src/utils.cpp ./src/config_param_parsing.cpp ./src/win_console_colors.cpp ./src/graph.cpp
set output=./bin/test.exe
set arguments=-Wall -pedantic -Wextra -Werror -O0 -g
set libraries=-lkernel32

set sanitize=%1
if %sanitize%==true set add_args=-fsanitize=address

clang %src% -o %output% %arguments% --include-directory src/include --include-directory src/include/scip %add_args% %libraries% -std=c++20