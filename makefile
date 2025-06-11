build: 
	cmake -Bbuild -H. -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -G Ninja

clean:
	rm -rf build
