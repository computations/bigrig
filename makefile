build: 
	cmake -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES

clean:
	rm -rf build bin tests/bin
