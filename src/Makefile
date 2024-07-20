CC = g++
CFLAGS = -c -Wall -Werror -Wextra -pedantic -std=c++17
LDFLAGS = -std=c++17
SOURCE = $(wildcard s21_*.cpp)
HEADER = $(wildcard s21_*.h)
OBJECTS = $(patsubst %.cpp,%.o, ${SOURCE})
TESTS_SOURCE = tests.cpp
TESTS_OBJECTS = $(patsubst %.cpp,%.o, ${TESTS_SOURCE})
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
  OPEN_CMD = xdg-open
	ADD_LIB = -lgtest -lsubunit -lm -lrt -lpthread -D_GNU_SOURCE
	BREW = .
endif
ifeq ($(UNAME_S),Darwin)
	OPEN_CMD = open
	ADD_LIB =
	BREW := ${HOME}/homebrew
endif

PATH := $(BREW)/bin:$(PATH)

all: s21_matrix_oop.a

s21_matrix_oop.a: ${SOURCE}
	$(CC) $(CFLAGS) $(ADD_LIB) $(SOURCE)
	ar rcs $@ $(OBJECTS)
	ranlib $@
	-rm -rf *.o

test: ${TESTS_SOURCE} s21_matrix_oop.a
	$(CC) $(LDFLAGS) -o test $^ -lgtest
	./test

gcov_report: ${SOURCE}
	$(CC) --coverage ${SOURCE} ${TESTS_SOURCE} -o gcov_test -lgtest -lgtest_main -std=c++17
	chmod +x gcov_test
	./gcov_test
	lcov -t "gcov_test" --ignore-errors inconsistent -o s21_test.info --no-external -c -d .
	genhtml --ignore-errors inconsistent -o report s21_test.info
	$(OPEN_CMD) ./report/index.html

lcov:
ifeq ("", "$(shell PATH=$(PATH) which lcov)")
	$(error Need to install lcov)
endif

valgrind: clean
	make s21_matrix_oop.a
	$(CC) -g ${TESTS_SOURCE} -std=c++17 s21_matrix_oop.a -o test -lgtest -lpthread -D_GNU_SOURCE 
	CK_FORK=no valgrind --trace-children=yes --track-fds=yes --track-origins=yes --leak-check=full --show-leak-kinds=all --verbose --log-file=RESULT_VALGRIND.txt ./test

style: clean
	cp ../materials/linters/.clang-format ../src/.clang-format
	clang-format -style=Google -n *.cpp *.h
	rm .clang-format

clean:
	-rm -rf *.o && rm -rf *.gcno
	-rm -rf *.a && rm -rf *.gcda
	-rm -rf *.info && rm -rf *.gcov
	-rm -rf ./gcov_test && rm -rf ./gcov_report
	-rm -rf ./report/
	-rm -rf ./test
	-rm -rf s21_test
	-rm -rf RESULT_VALGRIND.txt

rebuild: clean all
.PHONY: all clean check lcov