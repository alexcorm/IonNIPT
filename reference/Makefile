.PHONY: clean 


WISECONDOR_PATH=../wisecondor
RAW_PATH=raw
SRC= $(wildcard *.bam)
GCC = $(SRC:.bam=.gcc)
GC_COUNT=./ref/gccount
HG_19=/PUBLIC_DATA/ReferenceGenomes/UCSChg19/UCSC.hg19.fa 
IN_REF=./in/

all: reference
	 
%.pickle: %.bam
		  samtools rmdup -s $^ - | samtools view - -q 1 | python $(WISECONDOR_PATH)/consam.py -outfile $@

%.gcc: %.pickle $(GC_COUNT)
	   python $(WISECONDOR_PATH)/gcc.py $< $(GC_COUNT) $@


$(GC_COUNT): $(HG_19)
			 mkdir ./ref 
			 python $(WISECONDOR_PATH)/countgc.py $< $@

$(IN_REF): $(GCC)
		    mkdir $@
		    mv *.gcc $@

reference: $(IN_REF) 
			mkdir $(IN_REF)
		    python $(WISECONDOR_PATH)/newref.py $< $@

clean:
	@rm *.pickle
	@rm *.gcc 


