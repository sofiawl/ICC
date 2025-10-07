# Caminho base do LIKWID
LIKWID_HOME=/home/soft/likwid
LIKWID_INCLUDE=${LIKWID_HOME}/include
LIKWID_LIB=${LIKWID_HOME}/lib

# PROGRAMA
PROG=resolveEDO
OBJS=utils.o edo.o
VERIF=verificaEP02
CC=gcc

# Acrescentar onde apropriado as opções para incluir uso da biblioteca LIKWID 
CFLAGS=-O0 -DLIKWID_PERFMON -I${LIKWID_INCLUDE}
LFLAGS=-L${LIKWID_LIB} -llikwid -lm

# Lista de arquivos para distribuição
DISTFILES = *.c *.h LEIAME* Makefile
DISTDIR = ${USER}

.PHONY: all clean purge dist

%.o: %.c %.h
	$(CC) -c $(CFLAGS) -o $@ $<
all: $(PROG)

$(PROG): % : %.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean:
	@echo "Limpando sujeira ..."
	@rm -f *~ *.bak

purge:  clean
	@echo "Limpando tudo ..."
	@rm -f $(PROG) *.o core a.out $(DISTDIR) $(DISTDIR).tar

testeFormato: $(PROG) $(VERIF)
	@cat teste.dat | ./$(PROG) | ./$(VERIF)

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tar) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)

