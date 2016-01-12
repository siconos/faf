cmake_minimum_required(VERSION 2.8)

find_program(frama_c_command frama-c)

message(STATUS "frama-c command is: ${frama-c}")

execute_process(COMMAND ${frama_c_command} -print-path OUTPUT_VARIABLE frama_c_path)

function(frama_c_target name)

  set (file ${CMAKE_CURRENT_BINARY_DIR}/test_${name}.c)

  set(args -cpp-extra-args=-DFUNCODEGEN_CHECK -cpp-extra-args=-I${CMAKE_CURRENT_SOURCE_DIR} -cpp-extra-args=-I${frama_c_path}/libc -cpp-extra-args=-I${frama_c_path} -cpp-extra-args=-std=c99 -kernel-msg-key pp -val ${file} -print -ocode ${name}_cil.c -then -wp ${file} -no-print -wp-model nat+real  -wp-rte -wp-prover alt-ergo   -wp-out temp -wp-verbose 2  -then -report)


  add_custom_target(${name}_check
    COMMAND ${frama_c_command} ${args}
    DEPENDS ${name}_codegen
    COMMENT value analysis for ${name})

endfunction()

function(add_static_test)
  
  set(name ${ARGV0})
  set(exe ${ARGV1})
  set(args ${ARGV2})

  message(STATUS "adding test: ${name}")

  add_custom_target(${name}_codegen
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${exe}
    COMMENT generating C file: ${name}.c
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${exe} ${args} > ${CMAKE_CURRENT_BINARY_DIR}/test_${name}.c)


  frama_c_target(${name})

  add_test(${name} ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${name}_check)

  set_tests_properties(${name} PROPERTIES FAIL_REGULAR_EXPRESSION "Dead;Alarm;To be validated")

endfunction()

add_test(basic_generation ${CMAKE_CURRENT_SOURCE_DIR}/misc.py)
set_tests_properties(basic_generation PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR})
add_static_test(norm2d norm2d.py)
add_static_test(AlartCurnier ACJM.py 0)
add_static_test(JeanMoreau ACJM.py 1)

