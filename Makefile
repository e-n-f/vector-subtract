index: index.c
	cc -g -Wall -O3 -o index index.c jsonpull/jsonpull.c -lm

H = $(wildcard *.h) $(wildcard *.hpp)
C = $(wildcard *.c) $(wildcard *.cpp)

indent:
	clang-format -i -style="{BasedOnStyle: Google, IndentWidth: 8, UseTab: Always, AllowShortIfStatementsOnASingleLine: false, ColumnLimit: 0, ContinuationIndentWidth: 8, SpaceAfterCStyleCast: true, IndentCaseLabels: false, AllowShortBlocksOnASingleLine: false, AllowShortFunctionsOnASingleLine: false, SortIncludes: false}" $(C) $(H)
