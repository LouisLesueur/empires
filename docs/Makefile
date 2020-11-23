SOURCES = \
	index \

DEP = $(wildcard *.sty *.tex *.jpg *.png)
TARGETS = $(addprefix $(BUILDDIR),$(addsuffix .html,$(SOURCES)))

# Change LaTeX engine
PARAMETERS = --template GitHub.html5
PARAMETERS += --katex

all: $(TARGETS)

%.html : %.md $(DEP)
	pandoc $(PARAMETERS) $< -o $@

clean:
	@rm -f $(TARGETS)
