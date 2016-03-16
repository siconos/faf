cmake_minimum_required(VERSION 2.8)

include(CMakeParseArguments)
find_program(frama_c_command frama-c)

message(STATUS "frama-c command is: ${frama-c}")

execute_process(COMMAND ${frama_c_command} -print-path OUTPUT_VARIABLE frama_c_path)

function(frama_c_target name)

  set (file ${CMAKE_CURRENT_BINARY_DIR}/test_${name}.c)

  set(args -cpp-extra-args=-DFUNCODEGEN_CHECK -cpp-extra-args=-I${CMAKE_CURRENT_SOURCE_DIR} -cpp-extra-args=-I${frama_c_path}/libc -cpp-extra-args=-I${frama_c_path} -cpp-extra-args=-std=c99 -kernel-msg-key pp -val ${file} -slevel 16384 -print -ocode ${name}_cil.c -then -wp ${file} -no-print -wp-model nat+real  -wp-rte -wp-alt-ergo-opt="-backward-compat"  -wp-out temp -wp-verbose 2 -wp-prover alt-ergo  -then -report)


  add_custom_target(${name}_check
    COMMAND ${frama_c_command} ${args}
    DEPENDS ${name}_codegen
    COMMENT value analysis for ${name})

endfunction()

function(add_static_test)
  
  set(options)
  set(oneValueArgs NAME EXE)
  set(multiValueArgs ARGS)
  cmake_parse_arguments(add_static_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  
  set(name ${add_static_test_NAME})
  set(exe ${add_static_test_EXE})
  set(args ${add_static_test_ARGS})


  message(STATUS "adding ${name} test: ${exe} ${args}")

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
add_static_test(NAME norm2d EXE norm2d.py)
add_static_test(NAME AlartCurnier EXE ACJM.py ARGS AlartCurnier FAB)
add_static_test(NAME JeanMoreau EXE ACJM.py ARGS JeanMoreau FAB)

