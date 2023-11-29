CC=gcc
BUILDFLAGS=-std=c11 -Wall -Wextra -Werror
CFLAGS=$(BUILDFLAGS) -c

TESTFLAGS := -lcheck -lpthread -pthread

TESTFLAGS +=
ifeq ($(shell uname -s),Linux)
	TESTFLAGS += -lm -lrt -lsubunit
endif

GCOVFLAGS = -fprofile-arcs -ftest-coverage

TESTS_SOURCES=$(addprefix tests/, matrix_test.c)
SOURCES=s21_matrix.c s21_matrix_helpers.c
OBJECTS=$(SOURCES:.c=.o)

FILES_TO_CHECK=$(shell find . -name '*.c' && find . -name '*.h' )

MAINFILE=${MAIN}

.PHONY: test rebuild clean stylecheck fmt s21_matrix.a

all: clean s21_matrix.a

%.o: %.c
	$(CC) $(CFLAGS) -g $^ -o $@ 

s21_matrix.a: $(OBJECTS)
	ar rcs s21_matrix.a $^
	ranlib s21_matrix.a

gcov_report: test
	lcov -t test.out -o test.info -c -d . --no-external
	genhtml -o report test.info
	open report/index.html

gcovr_report: test
	mkdir gcovr/
	gcovr --html-details gcovr/report.details.html
	open gcovr/report.details.html

report_nohtml: test
	gcov *gcno > coverage.txt
	open coverage.txt

test: clean test_s21_matrix.a
	$(CC) -g $(BUILDFLAGS) $(TESTS_SOURCES) $(GCOVFLAGS) $(TESTFLAGS) -L. s21_matrix.a -o test.out
	./test.out

test_s21_matrix.a: all_functions
	ar rcs s21_matrix.a *.o
	ranlib s21_matrix.a

all_functions: $(SOURCES)
	$(CC) -g $(CFLAGS) $^ $(GCOVFLAGS)

rebuild: clean all

main: clean s21_matrix.a
	@[[ ! -z "$(MAINFILE)" ]] || (echo 'define MAIN env variable!' && exit 1)
	$(CC) -o main.out -g $(BUILDFLAGS) $(MAINFILE).c s21_matrix.a

clean: 
	@rm -rf  $(OBJECTS)
	@rm -rf report/ gcovr/ s21_matrix.a *.gcda *.gcno *.gcov *.out *.info *.o coverage.txt *.dSYM

.clang-format:
	@cp ../materials/linters/.clang-format .
	
stylecheck: .clang-format
	@clang-format -n $(FILES_TO_CHECK)
	
fmt: .clang-format
	@clang-format -i $(FILES_TO_CHECK)

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem $(FILES_TO_CHECK)
