@echo off
@REM -fsanitize=address -g3 -fno-omit-frame-pointer -fno-optimize-sibling-calls -target i386-pc-windows-msvc

SET src=./src/source.cpp ./src/utils.cpp ./src/config_param_parsing.cpp ./src/win_console_colors.cpp ./src/cycle_extraction.cpp ./src/ilp.cpp ./src/print/printing.cpp ./src/graph.cpp
SET output=./bin/ham.exe
SET arguments=-Wall -pedantic -Wextra -Ofast -g0
SET libraries=-lkernel32 -L./external/ -llibscip -lsoplex

@rem set sanitize=%1
@rem if %sanitize%==true set add_args=-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls

@REM cl %src% /O2 /I ./src/include /I ./src/include/scip /std:c++latest /EHsc /MD kernel32.lib /MD ./external/libscip32.lib /MD ./external/soplex32.lib

@REM del *.obj

@REM rename source.exe ilp_ls1.exe

@REM move ilp_ls1.exe bin/

clang %src% -o %output% %arguments% --include-directory src/include --include-directory src/include/scip %add_args% %libraries% -std=c++20