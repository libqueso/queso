all:	all_users all_developers

clean_all: clean_users clean_developers clean_exec

clean_users:
	cd users; make clean_all

clean_developers:
	cd developers; make clean_all

clean_exec:
	rm -f *~

all_users:
	cd users; make

all_developers:
	cd developers; make
