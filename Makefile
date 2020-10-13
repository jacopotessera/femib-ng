.PHONY: clean prepare cmake

clean:
	rm -rf build/

prepare:
	mkdir -p build

cmake: clean prepare
	cmake -H. -Bbuild

