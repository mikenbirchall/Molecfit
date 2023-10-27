GIT_URL=https://jugit.fz-juelich.de/mlz/libcerf.git
PKG_NAME=libcerf

what:
	@echo "all		-> Download libcerf package"
	@echo "libcerf		-> Download libcerf package"
	@echo "build_libcerf	-> Build downloaded libcerf package"
	@echo "clean_libcerf	-> Remove libcerf build"

help: what

all: $(PKG_NAME)


$(PKG_NAME):
	git clone $(GIT_URL)

clean_$(PKG_NAME):
	rm -rdf $(PKG_NAME)

build_$(PKG_NAME):$(PKG_NAME)
	cd $(PKG_NAME); mkdir build
	cd $(PKG_NAME)/build; cmake ..
	cd $(PKG_NAME)/build; make all
