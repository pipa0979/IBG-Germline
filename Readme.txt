valgrind --tool=memcheck program_name
valgrind --tool-name=memcheck --leak-check=yes program_name
valgrind --tool-name=memcheck --leak-check=yes --show-reachable=yes program_name