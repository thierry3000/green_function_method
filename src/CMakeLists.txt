add_executable(
green_function_method 
main.cpp
analyzenet.cpp
bicgstab.cpp
blood.cpp
contour.cpp
contr_lines.cpp
contr_shade.cpp
convect.cpp
eval.cpp
gaussj.cpp
greens.cpp
histogram.cpp
initgreens.cpp
input.cpp
ludcmp.cpp
nrutil.cpp
outboun.cpp
picturenetwork.cpp
putrank.cpp
readsources.cpp
setuparrays0.cpp
setuparrays1.cpp
setuparrays2.cpp
testconvect.cpp
tissrate.cpp
)
find_package(Boost REQUIRED COMPONENTS "system")
message("Boost_SYSTEM_LIBRARY: ${Boost_SYSTEM_LIBRARY}")
find_package(Boost REQUIRED COMPONENTS "filesystem")
message("Boost_FILESYSTEM_LIBRARY: ${Boost_FILESYSTEM_LIBRARY}")

target_link_libraries(green_function_method boost_filesystem boost_system)
