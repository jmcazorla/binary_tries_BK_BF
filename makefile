CFLAGS = -std=c++11 -O3 -march=native
CC = g++
OBJDIR = obj
BUILDIR = build
SDSLFLAGS = -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64
vpath %.cpp src
vpath %.hpp src

program: $(OBJDIR)/intersection_query_log.o \
         $(OBJDIR)/barbay_kenyon_log.o \
		 $(OBJDIR)/delta_barby_kenyon_log.o \
		 $(OBJDIR)/delta_trie_certificate_log.o
	$(CC) $(CFLAGS) test/build_tries.cpp -o $(BUILDIR)/build_tries.out $(SDSLFLAGS) 
	$(CC) $(CFLAGS) test/build_tries_h.cpp -o $(BUILDIR)/build.out $(SDSLFLAGS) 
	$(CC) -o $(BUILDIR)/queries_bk.out $(OBJDIR)/barbay_kenyon_log.o  $(CFLAGS) $(SDSLFLAGS)
	$(CC) -o $(BUILDIR)/queries_delta_bk.out $(OBJDIR)/delta_barby_kenyon_log.o  $(CFLAGS) $(SDSLFLAGS)
	$(CC) -o $(BUILDIR)/queries_delta_trie.out $(OBJDIR)/delta_trie_certificate_log.o  $(CFLAGS) $(SDSLFLAGS)
	$(CC) -o $(BUILDIR)/queries.out $(OBJDIR)/intersection_query_log.o  $(CFLAGS) $(SDSLFLAGS) -pthread 



$(OBJDIR)/barbay_kenyon_log.o: test/barbay_kenyon_log.cpp
	mkdir -p obj
	$(CC) -c -o $@ test/barbay_kenyon_log.cpp $(CFLAGS)  $(SDSLFLAGS)

$(OBJDIR)/delta_barby_kenyon_log.o: test/delta_barby_kenyon_log.cpp
	mkdir -p obj
	$(CC) -c -o $@ test/delta_barby_kenyon_log.cpp $(CFLAGS)  $(SDSLFLAGS)

$(OBJDIR)/delta_trie_certificate_log.o: test/delta_trie_certificate_log.cpp
	mkdir -p obj
	$(CC) -c -o $@ test/delta_trie_certificate_log.cpp $(CFLAGS)  $(SDSLFLAGS)

$(OBJDIR)/intersection_query_log.o: test/intersection_query_log.cpp
	mkdir -p obj
	$(CC) -c -o $@ test/intersection_query_log.cpp $(CFLAGS) $(SDSLFLAGS)

clean:
	rm -f core $(OBJDIR)/*.o
	rm -f core $(BUILDIR)/*.out