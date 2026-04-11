PKG    := TLPR
CXX    := /opt/homebrew/opt/llvm@17/bin/clang++
CXXFLAGS := -std=c++17 -O2 -DSTANDALONE_BUILD \
            -I inst/include

.PHONY: build install clear-cache lint fmt test-cpp

## Build the package shared library in-place (no install)
build:
	R CMD build .

## Install the package into the active R library
install:
	R CMD INSTALL --preclean .

## Remove all compiled artefacts and tooling caches
clear-cache:
	find src -maxdepth 1 \( -name "*.o" -o -name "*.so" -o -name "*.dll" -o -name "*.plist" \) -delete
	rm -rf src/.cache
	find . -maxdepth 1 -name "$(PKG)_*.tar.gz" -delete

## Build and run C++ unit tests (independent of R)
test-cpp:
	$(CXX) $(CXXFLAGS) src/hilbert.cpp tests/cpp/test_hilbert.cpp -o /tmp/tlpr_test_hilbert
	/tmp/tlpr_test_hilbert

## Run clang-tidy on all project C++ sources
lint:
	@for f in src/highs_model.cpp src/utils.cpp src/cartesian.cpp src/reader.cpp src/hilbert.cpp; do \
	  echo "--- $$f ---"; \
	  /opt/homebrew/opt/llvm@17/bin/clang-tidy -p src/ "$$f" 2>&1 \
	    | grep -E "^$$(pwd)/(src|inst)" || echo "clean"; \
	done

## Check formatting (non-destructive)
fmt:
	@/opt/homebrew/opt/llvm@17/bin/clang-format --dry-run --Werror \
	  src/highs_model.cpp src/utils.cpp src/cartesian.cpp src/reader.cpp src/hilbert.cpp \
	  inst/include/*.hpp \
	  && echo "Format OK" || echo "Format issues found — run: clang-format -i <file>"
