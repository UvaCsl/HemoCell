MAKE_PID := $(shell echo $$PPID)
JOB_FLAG := $(filter -j%, $(subst -j ,-j,$(shell ps T | grep "^\s*$(MAKE_PID).*$(MAKE)")))
JOBS     := $(subst -j,,$(JOB_FLAG))
ifeq ($(JOBS),)
	JOBS := "1"
endif


DIRS := $(shell find . -mindepth 1 -maxdepth 1 -type d)
CMAKES := $(DIRS:%=%/CMakeLists.txt)
EXECS := $(join $(DIRS), $(DIRS:./%=/%))

.PHONY: cmakefiles ../build/hemocell/libhemocell.a $(EXECS)
.NOTPARALLEL: ../build/hemocell/libhemocell.a

all: executables

executables: $(EXECS) 

cmakefiles: $(CMAKES)

$(CMAKES) : CMakeLists_template.txt 
	cp  $? $@
	sed -i 's/FOLDER_NAME__/$(@D)/g' $@


../build/hemocell/libhemocell.a:
	cd ../build/hemocell && cmake . && make -j $(JOBS)

$(EXECS) : ../build/hemocell/libhemocell.a cmakefiles
	cd $(@D) && mkdir -p build && cd build && cmake ../ && make -j1

clean:
	rm -rf $(CMAKES) $(EXECS) */build
