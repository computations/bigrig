extra_flags = -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -G Ninja

build:
	@cmake -Bbuild -H. $(extra_flags)

clang: extra_flags += -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang
clang: build

clean:
	rm -rf build
