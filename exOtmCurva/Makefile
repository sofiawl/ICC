# PROGRAMA
    PROGS = ajustePol gera_entrada ajustePolMEGAOTIM
    OBJS = utils.o

       LIKWID = /home/soft/likwid
 LIKWID_FLAGS = -I$(LIKWID)/include
  LIKWID_LIBS = -L$(LIKWID)/lib

    AVX_FLAGS = -mavx2 -march=native -fopt-info-vec

# Compilador
    CC     = gcc -O3 $(AVX_FLAGS) -Wno-unused-result 
    CFLAGS = -DLIKWID_PERFMON $(LIKWID_FLAGS) 
    LFLAGS = $(LIKWID_LIBS) -llikwid -lm

# Lista de arquivos para distribuição.
# LEMBRE-SE DE ACRESCENTAR OS ARQUIVOS ADICIONAIS SOLICITADOS NO ENUNCIADO DO TRABALHO
DISTFILES = *.c *.h LEIAME* Makefile 
DISTDIR = swl24

.PHONY: all clean purge dist

%o: %c utils.h
	$(CC) -o $@ $(CFLAGS) $^

all: $(PROGS)

$(PROGS) : % : %.o $(OBJS)
	$(CC) -o $@ $(CFLAGS) $^ $(LFLAGS)

clean:
	@echo "Limpando sujeira ..."
	@rm -f *~ *.bak core 

purge:  clean
	@echo "Limpando tudo ..."
	@rm -f $(PROGS) *.o a.out $(DISTDIR) $(DISTDIR).tar

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chvzf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
