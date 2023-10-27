GIT_URL=https://jugit.fz-juelich.de/mlz/libcerf.git
PKG_NAME=libcerf

all: $(PKG_NAME)


$(PKG_NAME):
	git clone $(GIT_URL)

clean_$(PKG_NAME):
	rm -rdf $(PKG_NAME)

build_$(PKG_NAME):$(PKG_NAME)
	cd $(PKG_NAME); mkdir build
	cd $(PKG_NAME)/build; cmake ..
	cd $(PKG_NAME)/build; make all
